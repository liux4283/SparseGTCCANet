# library(WGCNA)
# install.packages("igraph")
# install.packages('rlist')
# install.packages('rTensor')
# install.packages('purrr')
# install.packages('furrr')
# need to install GO.db, preprocessCore and impute package from bioconductor
#' All-at-once multi-omics network analysis pipeline design to implement SGTCCA-Net algorithm. This algorithm starts with 
#' running sparse generalized tensor canonical correlation analysis to extract the canonical weight, then construct adjacency matrices
#' based on canonical weights. After global adjacency matrices are constructed, the next step is to prune each adjacency matrix so that 
#' only important nodes are kept to construct multi-omics module, this is achieved by applying PageRank algorithm and evaluating through
#' correlation with respect to phenotype.
#'
#' @param data_list A list contain multi-omics data in matrix form, the rows of each data are subjects and the columns are 
#' molecular features, subject id need to be matched before creating the list.
#' @param correlation_list A list of correlation structure of interest for SGTCCA, for instance, if users are interested in 
#' the 4-way correlation between data 1-4, and the 3-way correlation between data 1, 2 and 4, then it has the form of list(c(1,2,3,4), c(1,2,4)).
#' @param num_solutions Number of solution need to extract for SGTCCA, equivlent to number of network modules.
#' @param prob_subsamp_prop Proportion of feature subsampled when calculating the covariance tensors/matrices density for SGTCCA subsampling.
#' @param prop_subsamp_num Total number of subsample used to estimate covariance tensors/matrices density.
#' @param common_subsamp_prop Proportion of of common subsampled feature for Turbo-SMT algorithm between iterations.
#' @param distinct_subsamp_prob Proportion of iteration-specific distinct subsampled features for each subsample iteration in Turbo-SMT.
#' @param num_workers Number of cores used for parallel computing, algorithm is parallizable with respect to number of subsamples.
#' @param num_iterations Total number of subsample iterations.
#' @param network_exclude_data: data (index) that are excluded from the network construction, usually this is the phenotype. For instance, if we want to exclude
#' the 4th dataset from the network construction, then network_exclude_data = 4.
#' @param data_type A vector specifying the type of data, and example would be c('gene','protein', 'metabolite').
#' @param summary_method Network summarization method used to summarize and trim network, either be PCA or NetSHy.
#' @param pheno A vector with phenotype data of interest to check the correlation between summarization score and phenotype for network trimming (currently only support single phenotype).
#' @param min_mod_size The minimum desirable module size after network trimming.
#' @param max_mod_size The maximum desirable module size after network trimming.
#' @param saving_dir directory where users want to save the result.


SGTCCA_Net <- function(data_list, correlation_list, num_solutions = 1, 
                       prob_subsamp_prop = 0.3, prob_subsamp_num = 200,
                       common_subsamp_prop = 0.08, distinct_subsamp_prop = 0.02,
                       num_workers = 5, num_iterations = 10, network_exclude_data = NULL,
                       data_type = NULL, summary_method = 'NetSHy', pheno = NULL,
                       min_mod_size = 20, max_mod_size = 300, saving_dir = NULL)
{
  # check if pheno argument is provided
  if (is.null(pheno))
    stop('phenotype data must be provided!')
  # check if saving directory is provided
  if (is.null(saving_dir))
    stop('saving directory must be provided!')
  # center, scale and transpose data
  data_list <- lapply(data_list, function(x){t(scale(x, center = TRUE, scale = TRUE))})
  
  # calculate sub-sampling density vector
  density <- subsamp_density(data_list = data_list, correlation_list = correlation_list, 
                             prob_subsamp_prop = prob_subsamp_prop, 
                             prob_subsamp_num = prob_subsamp_num, num_workers = num_workers
  ) 
  # print progress
  cat('###### Subsampling probability calculation completed! ######')
  
  # perform sgtcca algorithm
  gtcca_result <- turbo_smt(data_list, correlation_list, num_iterations, subsampling_density = density, 
                            common_subsamp_prop, distinct_subsamp_prop, 
                            num_workers, num_solutions)
  # print progress
  cat('###### SGTCCA completed! ######')
  
  # finding solution correspondence between subsamples
  feature_size <-  unlist(lapply(data_list, nrow))
  aggregated_network <- network_aggregation(gtcca_result, num_features = feature_size, num_solutions, data_length = length(data_list), correlation_length = length(correlation_list), correlation_list)
  
  
  
  # store aggregated weight (exclude dataset we don't wish to appear in the network module)
  cc_weight_network <- rlist::list.rbind(aggregated_network[[3]][-network_exclude_data])
  # extract feature labels and put it into the adjacency matrix.
  feature_label <- unlist(lapply(data_list, row.names))
  
  
  # network trimming
  types <- unlist(purrr::map((1:length(feature_size))[-network_exclude_data], function(x){
    rep(data_type[x], feature_size[x])
  }))
  
  # print progress
  cat('###### Network analysis completed! ######')
  
  # print(types)
  # combine data
  combined_data <- t(rlist::list.rbind(data_list[-network_exclude_data]))
  # calculate correlation matrix
  cor_mat <- cor(combined_data)
  # construct adjacency matrices
  abar <- list()
  for (j in 1:num_solutions)
  {
    abar[[j]] <- getAdj(cc_weight_network[,j], FeatureLabel = feature_label)
    # print(dim(abar[[j]]))
    # name global network
    names(abar[j]) <- paste0('Module ', j)
    # network trimming and summarization
    Network_Summarization_PPR_single(Abar = abar[[j]],CorrMatrix = cor_mat, type = types, data = combined_data, Pheno = pheno,
                                     ModuleIdx = j, min_mod_size = min_mod_size, max_mod_size = max_mod_size, method = summary_method, network_preference =  network_preference, saving_dir = saving_dir)
  }
  cat("############################")
  cat('SGTCCA-Net pipeline complete!')
  return(list(global_network = abar, correlation_matrix = cor_mat))
  
}  



# get adjacency matrix from canonical weight
getAdj <- function(Ws, FeatureLabel = NULL){
  
  
  if(is.null(dim(Ws))){
    Abar <- Matrix::Matrix(abs(Ws) %o% abs(Ws), sparse = TRUE)
  }else{
    b <- nrow(Ws)
    Abar <- matrix(0, nrow = b, ncol = b)
    for(ind in seq_len(ncol(Ws))){
      w <- abs(Ws[ , ind])
      A <- Matrix::Matrix(w %o% w, sparse = TRUE)
      Abar <- Abar + A
    }
  }
  
  diag(Abar) <- 0
  Abar <- Abar/max(Abar)
  
  
  colnames(Abar) <- rownames(Abar) <- FeatureLabel
  return(Abar)
}

cov_comp <- function(data_list, correlation_list)
{
  # require each data to be feature by subject
  # get number of features for each view
  num_features <- unlist(lapply(data_list, function(x){nrow(x)}))
  #data_list <- lapply(data_list, function(x){abs(x)})
  cov_tensor_list <- list()
  # compute covariance tensor/matrix
  for (i in 1:length(correlation_list))
  {
    
    # if even number of dimensions
    if (length(correlation_list[[i]]) %% 2 == 0)
    {  
      cov_tensor_list[[i]] <- array(0, dim = num_features[correlation_list[[i]]])
      for (j in 1:dim(data_list[[1]])[2])
      {
        # print(paste0('the ', j, 'th iteration'))
        C <- outer(data_list[[correlation_list[[i]][1]]][,j], data_list[[correlation_list[[i]][2]]][,j])
        if (length(correlation_list[[i]]) >= 3)
        {
          for (k in 3:length(correlation_list[[i]]))
          {
            C <- outer(C, data_list[[correlation_list[[i]][k]]][,j])
          }
        }  
        cov_tensor_list[[i]] <- cov_tensor_list[[i]] + C
      }
      cov_tensor_list[[i]] <- abs(cov_tensor_list[[i]]/dim(data_list[[1]])[2])
    }else{
      cov_tensor_list[[i]] <- array(0, dim = num_features[correlation_list[[i]]])
      for (k in 1:length(correlation_list[[i]]))
      {
        # define new data
        new_data_list <- data_list
        new_data_list[[correlation_list[[i]][k]]] <- abs(new_data_list[[correlation_list[[i]][k]]])
        temp_cov_tensor <- array(0, dim = num_features[correlation_list[[i]]])
        # print(new_data_list[[k]])
        for (j in 1:dim(data_list[[1]])[2])
        {
          
          
          C <- outer(new_data_list[[correlation_list[[i]][1]]][,j], new_data_list[[correlation_list[[i]][2]]][,j])
          if (length(correlation_list[[i]]) >= 3)
          {
            for (h in 3:length(correlation_list[[i]]))
            {
              C <- outer(C, new_data_list[[correlation_list[[i]][h]]][,j])
            }
          }  
          temp_cov_tensor <- temp_cov_tensor + C
          #cov_tensor_list[[i]] <- cov_tensor_list[[i]] + C
        }
        cov_tensor_list[[i]] <- cov_tensor_list[[i]] + abs(temp_cov_tensor)
      }  
      cov_tensor_list[[i]] <- cov_tensor_list[[i]]/(length(correlation_list[[i]]) * dim(data_list[[1]])[2])
      #print(length(correlation_list[[i]]) * dim(data_list[[1]])[2])
    }
    
  }
  return(cov_tensor_list)
}

subsamp_density <- function(data_list,correlation_list, prob_subsamp_prop, prob_subsamp_num, num_workers)
{
  # making clusters for multisession computation
  future::plan(future::multisession, workers = num_workers)
  
  # execute multisession computation
  subsamp_clusters <- furrr::future_map(1:prob_subsamp_num, function(xx) {
    library(parallel)
    library(magrittr)
    library(future)
    # obtain subsampling index from each data type 
    subsamp_index <- purrr::map(data_list, function(x){sort(sample(1:nrow(x), ceiling(nrow(x) * prob_subsamp_prop), 
                                                                   replace = FALSE))})
    # subsampling indexing for each dataset and obtain subsample data
    subsamp_data <- purrr::map(1:length(data_list), function(x){
      as.matrix(data_list[[x]][subsamp_index[[x]],])
    })
    
    # in case there are dataset with one row
    subsamp_data <- lapply(subsamp_data, function(x){if(ncol(x) == 1) {t(x)}else{x} 
    })
    
    # obtain covariance component based on correlation of interest
    cov_sub <- cov_comp(subsamp_data, correlation_list)
    cov_sub <- lapply(cov_sub, abs)
    
    
    
    # initialize probability sampling vector
    prob_samp_vec <- purrr::map(1:length(data_list), function(x) 
      matrix(0, nrow = nrow(data_list[[x]]), ncol = 1))
    
    # initialize all-one vector for density calculation
    num_features <- unlist(lapply(subsamp_data, nrow))
    one <- purrr::map(num_features, function(x){
      t(rep(1,x))
    })
    
    # in case there are view with one feature
    one <- lapply(one, function(x){as.matrix(x)
    })
    
    for (i in 1:length(correlation_list))
    {
      # store all mode index
      all_mode <- 1:length(correlation_list[[i]])
      # store view index
      view_index <- correlation_list[[i]]
      # calculate probability vector
      for (j in 1:length(correlation_list[[i]]))
      {
        # temporarily store the unscale probability vector
        #temp <- as.vector(rTensor::modeSum(rTensor::as.tensor(cov_sub[[i]]), 
        #all_mode[!all_mode %in% j], drop = FALSE)@data)
        temp <- as.vector(rTensor::ttl(rTensor::as.tensor(cov_sub[[i]]), one[view_index[!view_index %in% view_index[j]]],all_mode[!all_mode %in% j])@data)
        # add probability vector to corresponding view
        prob_samp_vec[[view_index[j]]][subsamp_index[[view_index[[j]]]],1] <- prob_samp_vec[[view_index[j]]][subsamp_index[[view_index[[j]]]],1] + 
          temp^2/sum(temp^2)
      }
      
    }
    gc()
    return(prob_samp_vec)
    
  },.progress = TRUE,.options = furrr::furrr_options(seed = TRUE))
  
  # convert nested list to list of data frames
  subsamp_prob_matrix_list <- purrr::map(1:length(data_list), 
                                         function(x){do.call(cbind.data.frame, 
                                                             lapply(subsamp_clusters, '[[', x))})
  
  
  # aggregate each subsampling probability matrix
  subsamp_prob <- purrr::map(1:length(subsamp_prob_matrix_list), function(x) 
  {
    apply(subsamp_prob_matrix_list[[x]], 2, function(y){y/pracma::Norm(y)})
    apply(subsamp_prob_matrix_list[[x]], 1, function(y){mean(y[y!=0])})
    
  })
  return(subsamp_prob)
}



# Writing SGTCCA algorithm (extract 1 solution only)
gtcca_update <- function(sub_data_list, correlation_list, alpha = 10000, deflation = TRUE, cov_residual = NULL, double_sparse = FALSE, beta = 2)
{
  # calculate covariance tensor/matrices (check if this is the first solution or not)
  if (deflation == TRUE)
  {
    cov_tensor_list <- cov_comp(data_list = sub_data_list, correlation_list = correlation_list)
    # convert to tensor
    cov_tensor_list <- lapply(cov_tensor_list, function(x){rTensor::as.tensor(x)})
  }else{
    # if that's for solution other than the first one, then use the covariance tensor residual instead
    cov_tensor_list <- cov_residual
  }
  
  
  # get each chunk length from the canonical weight vectors
  chunk_length_cc_weight <- unlist(lapply(sub_data_list, function(x){nrow(x)}))
  # get each chunk length from covariance parameters
  chunk_length_cov_par <- rep(1, length(correlation_list))
  # combine chunk length
  chunk_length <- c(chunk_length_cc_weight, chunk_length_cov_par)
  
  # source code of superdiagonal tensor
  superdiagonal_tensor <- function(num_modes,len,elements=1L){
    modes <- rep(len,num_modes)
    arr <- array(0, dim = modes)
    if(length(elements)==1) elements <- rep(elements,len)
    for (i in 1:len){
      txt <- paste("arr[",paste(rep("i", num_modes),collapse=","),"] <- ", elements[i],sep="")
      eval(parse(text=txt))
    }
    rTensor::as.tensor(arr)
  }
  
  # update rule (gradient)
  
  # optimization function
  opt_fun <- function(vec)
  {
    # partition by group (create group)
    chunk_group <- unlist(purrr::map(1:length(chunk_length), function(x){
      rep(x, chunk_length[x])
    }))
    # split long vector into different chunks
    par_list <- split(vec, chunk_group)
    # reconstruct rank-1 tensors
    cov_recon_list <- purrr::map(1:length(chunk_length_cov_par), function(x){
      sup_lambda <- superdiagonal_tensor(cov_tensor_list[[x]]@num_modes, 1, 
                                         par_list[[length(chunk_length_cc_weight)+x]])
      # extract a list of canonical weight vectors based on correlation_list
      par_cc_weight <- par_list[correlation_list[[x]]]
      # convert them into matrix for tensor-based calculation
      par_cc_weight <- lapply(par_cc_weight, function(x){matrix(x, ncol = 1)})
      # reconstruct rank-1 tensor
      recon <-  rTensor::ttl(sup_lambda,par_cc_weight,ms=1:cov_tensor_list[[x]]@num_modes)
      return(recon)
    })
    # get norm of each canonical weight vector
    norm_vec <- unlist(purrr::map(1:length(chunk_length_cc_weight), function(x){
      pracma::Norm(par_list[[x]])
    }))
    
    # optimization function construction (error + norm)
    opt_error <- sum(unlist(purrr::map(1:length(cov_tensor_list), function(x){
      rTensor::fnorm(cov_tensor_list[[x]] - cov_recon_list[[x]])^2
    })))
    opt_norm <- alpha * sum((norm_vec-1)^2)
    
    if (double_sparse == TRUE)
    {
      sparse_pen_loss <- unlist(purrr::map(1:length(chunk_length_cc_weight), function(x){
        sum(sqrt(par_list[[x]]^2 + 0.0001))
      }))
      opt_norm <- opt_norm + beta * sum(sparse_pen_loss)
    }
    
    # return overall function construction
    return(opt_error + opt_norm)
  }  
  
  # derivative function
  der_fun <- function(vec)
  {
    # partition by group (create group)
    chunk_group <- unlist(purrr::map(1:length(chunk_length), function(x){
      rep(x, chunk_length[x])
    }))
    # split long vector into different chunks
    par_list <- split(vec, chunk_group)
    # reconstruct rank-1 tensors
    cov_recon_list <- purrr::map(1:length(chunk_length_cov_par), function(x){
      sup_lambda <- superdiagonal_tensor(cov_tensor_list[[x]]@num_modes, 1, 
                                         par_list[[length(chunk_length_cc_weight)+x]])
      # extract a list of canonical weight vectors based on correlation_list
      par_cc_weight <- par_list[correlation_list[[x]]]
      # convert them into matrix for tensor-based calculation
      par_cc_weight <- lapply(par_cc_weight, function(x){matrix(x, ncol = 1)})
      # reconstruct rank-1 tensor
      recon <-  rTensor::ttl(sup_lambda,par_cc_weight,ms=1:cov_tensor_list[[x]]@num_modes)
      return(recon)
    })
    # get norm of each canonical weight vector
    norm_vec <- unlist(purrr::map(1:length(chunk_length_cc_weight), function(x){
      pracma::Norm(par_list[[x]])
    }))
    
    # calculate derivative for canonical weight vector
    der_cw <- purrr::map(1:length(sub_data_list), function(x)
    {
      # check which covariance component has a specific data type
      data_type_ind <- unlist(lapply(correlation_list, function(y){x %in% y}))
      # subset only those correlation structure
      sub_correlation_list <- correlation_list[data_type_ind]
      # check the mode index of a specific data type
      sub_correlation_index <- unlist(lapply(sub_correlation_list, function(z){
        which(z %in% x)
      }))
      # select only the subset of covariance tensor 
      cov_subset <- cov_tensor_list[data_type_ind]
      # select only the subset of reconstructed covariance tensor
      cov_recon_subset <- cov_recon_list[data_type_ind]
      # calculate gradient for canonical weight
      gradient_cw <- rowSums(do.call(cbind.data.frame,
                                     purrr::map(1:length(sub_correlation_index), function(z){
                                       # create a list of canonical weight for gradient calculation, this should be all canonical weight vector that are not taken gradient at this stage
                                       gradient_index <- sort(sub_correlation_list[[z]][!sub_correlation_list[[z]] %in% x], decreasing = TRUE)
                                       # canonical weight matrix list
                                       cw_list <- purrr::map(gradient_index, function(zz){
                                         as.matrix(par_list[[zz]])
                                       })
                                       if (length(cw_list) > 1)
                                       {(rTensor::k_unfold(cov_recon_subset[[z]], sub_correlation_index[z])@data - rTensor::k_unfold(cov_subset[[z]], sub_correlation_index[z])@data) %*% (par_list[[length(sub_data_list) + z]] * rTensor::khatri_rao_list(cw_list))} else{
                                         (rTensor::k_unfold(cov_recon_subset[[z]], sub_correlation_index[z])@data - rTensor::k_unfold(cov_subset[[z]], sub_correlation_index[z])@data) %*% (par_list[[length(sub_data_list) + z]] * cw_list[[1]])
                                       }
                                     })
      ))
      return(gradient_cw)
      
    })
    # plus shrinkage towards norm of 1
    der_cw <- purrr::map(1:length(sub_data_list), function(x){
      der_cw[[x]] + alpha * (par_list[[x]] - par_list[[x]]/pracma::Norm(par_list[[x]]))
    })
    
    # double sparse or not 
    if (double_sparse == TRUE)
    {
      der_cw <- purrr::map(1:length(sub_data_list), function(x){
        der_cw[[x]] + (beta/2) * (par_list[[x]]/sqrt(par_list[[x]]^2 + 0.0001))
      })
    }
    
    # calculate derivative for scaling parameters
    der_sp <- purrr::map(1:length(correlation_list), function(x){
      # extract a list of canonical weight vectors based on correlation_list
      par_cc_weight <- par_list[correlation_list[[x]]]
      # convert them into matrix for tensor-based calculation
      par_cc_weight <- lapply(par_cc_weight, function(z){matrix(z, nrow = 1)})
      # calculate gradient
      gradient_sp <- rTensor::ttl(cov_recon_list[[x]] - cov_tensor_list[[x]], par_cc_weight, 1:length(correlation_list[[x]]))@data
      return(gradient_sp)
    })
    
    # combine gradient into one single vector 
    gradient <- c(unlist(der_cw), unlist(der_sp))
    return(gradient)
  }
  
  # calculate the total length of the gradient vector
  length_vector <- sum(unlist(lapply(sub_data_list, nrow))) + length(correlation_list)
  # randomly initialize parameter
  vec <- rnorm(length_vector, 0, 1)
  
  # optimize objective function
  result <- Rcgmin::Rcgmin(par = vec, fn = opt_fun, gr = der_fun, control = list(maxit = 500))
  
  # partition by group (create group)
  chunk_group <- unlist(purrr::map(1:length(chunk_length), function(x){
    rep(x, chunk_length[x])
  }))
  # split long vector into different chunks
  cc_result <- split(result[['par']], chunk_group)
  
  # calculate reconstructed covariance tensor
  cov_recon_list <- purrr::map(1:length(chunk_length_cov_par), function(x){
    sup_lambda <- superdiagonal_tensor(cov_tensor_list[[x]]@num_modes, 1, 
                                       cc_result[[length(chunk_length_cc_weight)+x]])
    # extract a list of canonical weight vectors based on correlation_list
    par_cc_weight <- cc_result[correlation_list[[x]]]
    # convert them into matrix for tensor-based calculation
    par_cc_weight <- lapply(par_cc_weight, function(x){matrix(x, ncol = 1)})
    # reconstruct rank-1 tensor
    recon <-  rTensor::ttl(sup_lambda,par_cc_weight,ms=1:cov_tensor_list[[x]]@num_modes)
    return(recon)
  })
  # calculate residual covariance tensor
  cov_residual_list <- purrr::map(1:length(chunk_length_cov_par), function(x){
    cov_tensor_list[[x]] - cov_recon_list[[x]]
  })
  # calculate original covariance tensor for norm explained calculation
  if (deflation == FALSE)
  {
    # calculate original covariance tensor (overwrite to reduce memory cost)
    cov_tensor_list <- cov_comp(data_list = sub_data_list, correlation_list = correlation_list)
    # convert to tensor
    cov_tensor_list <- lapply(cov_tensor_list, function(x){rTensor::as.tensor(x)})
  }  
  # calculate percentage of norm explained
  per_norm_explained <- purrr::map(1:length(chunk_length_cov_par), function(x){
    1-rTensor::fnorm(cov_residual_list[[x]])/rTensor::fnorm(cov_tensor_list[[x]])
  })
  
  return(list(cc_result, cov_residual_list, per_norm_explained, der_fun(result[['par']]), result[['counts']]))
  
}


# Deflation Method 
gtcca_deflation <- function(data_list, correlation_list, num_solution)
{
  # source code of superdiagonal tensor
  superdiagonal_tensor <- function(num_modes,len,elements=1L){
    modes <- rep(len,num_modes)
    arr <- array(0, dim = modes)
    if(length(elements)==1) elements <- rep(elements,len)
    for (i in 1:len){
      txt <- paste("arr[",paste(rep("i", num_modes),collapse=","),"] <- ", elements[i],sep="")
      eval(parse(text=txt))
    }
    rTensor::as.tensor(arr)
  }
  # initialize canonical weight and parameter matrices
  num_features <- lapply(data_list, nrow)
  cc_weight_matrix <- purrr::map(1:length(data_list), function(x){
    matrix(0, nrow = num_features[[x]], ncol = num_solution)
  })
  
  # initialize canonical weight scaling number
  cc_scale_matrix <- matrix(0, nrow = length(correlation_list), ncol  = num_solution)
  
  # initialize norm percentage matrix
  norm_explained <- cum_norm_explained <- matrix(0, nrow = length(correlation_list), ncol = num_solution)
  # deflation starts
  i <- 1
  while (i  <= num_solution)
  {
    # run GTCCA
    if (i == 1)
    {
      cc_result_temp <- gtcca_update(data_list, correlation_list, deflation = TRUE)
    }else{
      cc_result_temp <- gtcca_update(data_list, correlation_list, deflation = FALSE, cov_residual = cc_result[[2]])
    }  
    
    # store result
    for (j in 1:length(data_list))
    {
      cc_weight_matrix[[j]][,i] <- cc_result_temp[[1]][[j]]
    }
    # unlist all scaling parameters
    scale_par <- unlist(cc_result_temp[[1]][(length(data_list)+1):(length(data_list) + length(correlation_list))])
    # store cumulative norm explained 
    cum_norm_explained[,i] <- unlist(cc_result_temp[[3]])
    # print(scale_par)
    # store result
    cc_scale_matrix[,i] <- scale_par
    
    # check if algorithm fails at certain iteration, if fails, rerun
    if ((cc_result_temp[[5]][2] < 20) || (sum(cum_norm_explained[,i]) < 0.01))
    {
      i <- i - 1
    }else{
      cc_result <- cc_result_temp
    }
    i <- i + 1
    
  }
  # calculate norm explained from percentage of cumulative norm explained
  norm_explained[,1] <- cum_norm_explained[,1]
  if (num_solution > 1)
  {  
    for (j in 2:ncol(cum_norm_explained))
    {
      norm_explained[,j] <- cum_norm_explained[,j] - cum_norm_explained[,j-1]
    }
  }  
  return(list(canonical_weight = cc_weight_matrix, scaling_parameters = cc_scale_matrix, residual_covariance_tensor = cc_result[[2]], norm_explained))
  
  
}

# turbo-smt algorithm
turbo_smt <- function(data, correlation_list, num_iterations, subsampling_density, common_subsamp_prop, distinct_subsamp_prop, num_workers, num_solutions, adjacency_measure = 'adjmatrix')
{
  # initialize common subsample index
  #common_index <- purrr::map(1:length(data), function(x){
  #  sample(1:nrow(data[[x]]), ceiling(nrow(data[[x]]) * common_subsamp_prop),
  #         replace = FALSE, prob = subsampling_density[[x]])
  #})
  common_index <- purrr::map(1:length(data), function(x){
    order(subsampling_density[[x]], decreasing = TRUE)[1:(ceiling(nrow(data[[x]]) * common_subsamp_prop))]
  })
  
  # making clusters for multisession computation
  future::plan(future::multisession, workers = num_workers)
  # start the turbo algorithm
  subsamp_clusters <- furrr::future_map(1:num_iterations, function(xx) {
    # distinct subsample index
    
    # index except for common 
    distinct_index <- purrr::map(1:length(data), function(x){
      if (nrow(data[[x]]) > 1)
      {
        index_left = (1:nrow(data[[x]]))[-common_index[[x]]]
        return(sample(index_left, ceiling(nrow(data[[x]]) * distinct_subsamp_prop),
                      replace = FALSE, prob = subsampling_density[[x]][-common_index[[x]]]))
      }else{
        return(sample(1:nrow(data[[x]]), ceiling(nrow(data[[x]]) * distinct_subsamp_prop),
                      replace = FALSE, prob = subsampling_density[[x]]))
      }
      
    })
    # merge two set of index
    all_index <- purrr::map(1:length(data), function(x){
      sort(unique(c(common_index[[x]], distinct_index[[x]])))
    })
    # sample the subset of data
    subset_data <- purrr::map(1:length(data), function(x){
      as.matrix(data[[x]][all_index[[x]],])
    })
    
    # in case there are dataset with one row
    subset_data <- lapply(subset_data, function(x){if(ncol(x) == 1) {t(x)}else{x} 
    })
    
    # run gtcca
    result <- try(gtcca_deflation(data_list = subset_data, correlation_list = correlation_list, num_solution = num_solutions))
    gc()
    return(c(result, sub_sample_index = list(all_index), shared_index = list(common_index)))
    
  },.progress = TRUE,.options = furrr::furrr_options(seed = TRUE))
  return(subsamp_clusters)
}




# aggregate canonical weight into different network modules 
network_aggregation <- function(turbo_result, num_features, num_solutions, data_length, correlation_length, correlation_list, is.trait = FALSE)
{
  turbo_result <- turbo_result[sample(1:length(turbo_result))]
  # filter out unsuccessful run
  successful_index <- unlist(lapply(turbo_result, length))
  turbo_result <- turbo_result[successful_index > 3]
  # check norm explained
  norm_check <- rep(0, length(turbo_result))
  for (i in 1:length(turbo_result))
  {
    per_norm_explained_sum <- apply(turbo_result[[i]][[4]],2,sum)
    if ((max(per_norm_explained_sum) < 0.01) || (max(per_norm_explained_sum) > 10))
    {
      norm_check[i] <- FALSE
    }else{
      norm_check[i] <- TRUE
    }
    
  }
  
  turbo_result <- turbo_result[which(norm_check == 1)]
  
  # create empty matrix for storing canonical weight
  cc_weight <- purrr::map(num_features, function(x){
    array(0, dim = c(x, num_solutions, length(turbo_result)))
  })
  # create empty matrix for storing percentage norm explained and scaling parameter
  per_norm_explained <- array(0, dim = c(correlation_length, num_solutions, length(turbo_result)))
  scaling_parameter <- array(0, dim = c(correlation_length, num_solutions, length(turbo_result)))
  
  
  # stitching the sub-sample canonical weight
  for (i in 1:length(turbo_result))
  {
    # store percentage norm explained
    if (!is.null(nrow(turbo_result[[i]][[4]])))
    {  
      per_norm_explained_sum <- apply(turbo_result[[i]][[4]],2,sum)
      turbo_result[[i]][[4]] <- turbo_result[[i]][[4]][,order(per_norm_explained_sum, decreasing = TRUE)]
      turbo_result[[i]][["scaling_parameters"]] <- turbo_result[[i]][["scaling_parameters"]][,order(per_norm_explained_sum, decreasing = TRUE)]
    }else{
      per_norm_explained_sum <- sum(turbo_result[[i]][[4]])
      turbo_result[[i]][[4]] <- turbo_result[[i]][[4]][order(per_norm_explained_sum, decreasing = TRUE)]
      turbo_result[[i]][["scaling_parameters"]] <- turbo_result[[i]][["scaling_parameters"]][order(per_norm_explained_sum, decreasing = TRUE)]
    }  
    
    per_norm_explained[,,i] <- turbo_result[[i]][[4]]
    scaling_parameter[,,i] <- turbo_result[[i]][["scaling_parameters"]]
    # store canonical weight
    for (j in 1:data_length)
    {
      if (num_solutions != 1)
      {
        turbo_result[[i]][["canonical_weight"]][[j]] <- as.matrix(turbo_result[[i]][["canonical_weight"]][[j]][,order(per_norm_explained_sum, decreasing = TRUE)])
      }
      cc_weight[[j]][turbo_result[[i]][["sub_sample_index"]][[j]],,i] <- abs(turbo_result[[i]][["canonical_weight"]][[j]])
    }
  }
  # finding column correspondence
  if (num_solutions != 1)
  {  
    for (i in 2:length(turbo_result))
    {
      # calculate correlation matrix between iterations based on shared sub-sample index
      temp_correlation <- abs(cor(cc_weight[[1]][turbo_result[[1]][['shared_index']][[1]], , 1], cc_weight[[1]][turbo_result[[i]][['shared_index']][[1]], , i]))
      #print(temp_correlation)
      # create empty vector to store solution correspondence
      temp_solution_correspondence <- c()
      # loop over to find the solution correspondence
      for (j in 1:num_solutions)
      {
        max_index <-  which.max(temp_correlation[j, ])
        
        temp_solution_correspondence <- c(temp_solution_correspondence, max_index)
        temp_correlation[, max_index] <- 0
      }
      # change the order of solution based on solution correspondence
      for (k in 1:data_length)
      {
        cc_weight[[k]][,,i] <- cc_weight[[k]][,temp_solution_correspondence,i]
      }
      # change the order of norm percentage explained
      per_norm_explained[,,i] <- per_norm_explained[,temp_solution_correspondence,i]
      # print(temp_solution_correspondence)
    }
  }  
  # normalize each canonical weight vector
  for (i in 1:length(turbo_result))
  {
    for(k in 1:data_length)
    {
      cc_weight[[k]][,,i] <- apply(as.matrix(cc_weight[[k]][,,i]), 2, function(x){x/pracma::Norm(x[turbo_result[[1]][["shared_index"]][[k]]])})
    }
  }
  
  # weight averaging
  # initialize weight averaged matrix
  cc_matrix <- purrr::map(num_features, function(x){
    array(0, dim = c(x, num_solutions))
  })
  
  # aggregate for canonical weight and scaling parameters
  for (i in 1:num_solutions)
  {
    for (j in 1:data_length)
    {
      #correlation_filter <- which(cor(cc_weight[[j]][,i,],cc_weight[[j]][,i,1]) > 0.8)
      #print(correlation_filter)
      if (num_features[j] != 1)
      {  
        cc_matrix[[j]][,i] <- apply(as.matrix(cc_weight[[j]][,i,]), 1, function(x){ifelse(sum(x) > 0, mean(x[x!=0]), 0)}) 
      }else{
        cc_matrix[[j]][,i] <- ifelse(sum(cc_weight[[j]][,i,]) > 0, mean(sum(cc_weight[[j]][,i,])[sum(cc_weight[[j]][,i,])!=0]), 0)
      }  
    }
  }
  
  norm_explained <- matrix(0, nrow = correlation_length, ncol = num_solutions)
  # scaling parameters 
  for (i in 1:num_solutions)
  {
    if (correlation_length != 1)
    {  
      norm_explained[,i] <- apply(per_norm_explained[,i,], 1, mean)
    }else{
      norm_explained[,i] <- mean(per_norm_explained[1,i,])
    }
  }
  
  
  cc_matrix <- lapply(cc_matrix, function(x){
    apply(x, 2, function(y){
      y/pracma::Norm(y)
    })
  })
  
  
  
  # aggregate networks '
  
  # initialize affinity matrix
  if (is.trait)
  { 
    # initialize affinity matrix
    affinity_matrix <- matrix(0, nrow = sum(num_features), ncol = sum(num_features))
    # get all possible pairwise structure
    combination <- combn(data_length, 2)
    # get cumulative sum for number of features 
    cum_ind <- c(0,cumsum(num_features[1:data_length]))
    for (i in 1:ncol(combination))
    {
      # check whether specific combination is in a certain correlation structure
      check_structure <- unlist(lapply(correlation_list, function(x){
        all(combination[,i] %in% x)
      }))
      
      for (j in 1:num_solutions)
      {
        # obtain the maximum scaling factor
        max_scaling <- max(norm_explained[check_structure,j]) 
        # obtain corresponding affinity matrix index
        affinity_matrix[(cum_ind[combination[1,i]] + 1):cum_ind[combination[1,i]+1], 
                        (cum_ind[combination[2,i]] + 1):cum_ind[combination[2,i]+1]] <- max_scaling * outer(cc_matrix[[combination[1,i]]][,j],cc_matrix[[combination[2,i]]][,j])  
      }  
      
      
    }
    
  }else{
    # reduce data length
    reduced_data_length <- data_length - 1
    
    # initialize affinity matrix
    affinity_matrix <- matrix(0, nrow = sum(num_features[1:reduced_data_length]), ncol = sum(num_features[1:reduced_data_length]))
    # get all possible pairwise structure
    combination <- combn(reduced_data_length, 2)
    # get cumulative sum for number of features 
    cum_ind <- c(0,cumsum(num_features[1:reduced_data_length]))
    
    # only retain the best canonical weight
    #cc_matrix <- lapply(cc_matrix, function(x){
    #  t(apply(x, 1, function(y){
    #    ifelse(y == max(y), y, 0)
    #  }))
    #})
    # print(dim(cc_matrix[[1]]))
    for (i in 1:ncol(combination))
    {
      # check whether specific combination is in a certain correlation structure
      check_structure <- unlist(lapply(correlation_list, function(x){
        all(combination[,i] %in% x)
      }))
      
      for (j in 1:num_solutions)
      {
        # obtain the maximum scaling factor
        max_scaling <- max(norm_explained[check_structure,j]) 
        # obtain corresponding affinity matrix index
        affinity_matrix[(cum_ind[combination[1,i]] + 1):cum_ind[combination[1,i]+1], 
                        (cum_ind[combination[2,i]] + 1):cum_ind[combination[2,i]+1]] <-affinity_matrix[(cum_ind[combination[1,i]] + 1):cum_ind[combination[1,i]+1], 
                                                                                                       (cum_ind[combination[2,i]] + 1):cum_ind[combination[2,i]+1]] + 
          max_scaling * outer(cc_matrix[[combination[1,i]]][,j],cc_matrix[[combination[2,i]]][,j])  
      }  
      
      
    }
  }
  
  
  
  return(list(cc_weight, per_norm_explained, cc_matrix, norm_explained, scaling_parameter, affinity_matrix))
}





selection <- function(corr_vec,corr_pheno,default_size = 10, cor_cut = 0.8,
                      network_preference = 'small', tol = 0.9)
{
  # Step 1: find all candidate within the customized range of the maximum correlation with respect to phenotype
  # Step 2: find the maximum/minimum size with cor at least cor_cut
  
  candidate_size_1 <- which(corr_pheno > (tol * max(corr_pheno)))
  candidate_size_2 <- which(corr_vec >= cor_cut)
  
  # find the maximum of the intersection
  if (network_preference == 'small')
    candidate <- min(intersect(candidate_size_1, candidate_size_2))
  if (network_preference == 'large')
    candidate <- max(intersect(candidate_size_1, candidate_size_2))
  
  # return the result 
  return (candidate + default_size - 1)
  
}



hybrid_score = function(X, A, is_alpha = TRUE, npc=1){
  #X: data matrix (n, p)
  #A: corresponding adjacency matrix
  #pc_id: PC index
  g = igraph::graph_from_adjacency_matrix(A, mode = "undirected", weighted = TRUE, diag = FALSE)
  # laplacian
  L2 <- igraph::graph_from_adjacency_matrix(A, mode = "undirected", weighted = TRUE,diag = FALSE) 
  L2 <- as.matrix(igraph::graph.laplacian(L2, normalized = F))
  # TOM
  TOM_A = WGCNA::TOMsimilarity(as.matrix(A), verbose = FALSE)
  
  alpha = igraph::graph.density(g)
  X= scale(X, center = TRUE, scale = TRUE)
  # weighted approach
  if (is_alpha == TRUE)
    temp  = (1-alpha)*(X %*% L2) + alpha*(X %*% TOM_A)
  else
    temp = (X %*% L2)
  temp = summary(prcomp(temp))
  h_score = temp$x[,1:npc] 
  importance = temp$importance[,1:npc]
  loading = temp$rotation[,1:npc]
  return(list(h_score, importance, loading))
  
}

Network_Summarization_PPR_single <- function(Abar, CorrMatrix, data,
                                             Pheno, type, ModuleIdx,
                                             min_mod_size, max_mod_size, damping = 0.9,
                                             method = 'NetSHy',network_preference, saving_dir){
  
  
  
  P1 <- p <- ncol(Abar)
  # Trim module by PPR
  net_ppr <- igraph::graph_from_adjacency_matrix(Abar, weighted = TRUE,
                                                 diag = FALSE, mode = "undirected")
  igraph::set_vertex_attr(net_ppr, "type", index = igraph::V(net_ppr), as.factor(type))
  # All parameter setups are based on the recommendation
  ranking <- igraph::page_rank(net_ppr, directed = FALSE, damping = damping, 
                               options = list(niter = 10^5,eps = 1e-06))
  # Obtain ranked protein names
  rank_names <- names(sort(ranking$vector, decreasing = TRUE))
  rank_value <- sort(ranking$vector, decreasing = TRUE)
  # Choose to include only the top proteins as we desired
  summary_correlation <- c()
  for (i in min_mod_size : max_mod_size)
  { 
    # print(paste0('iteration:', i))
    newM.node <- which(colnames(Abar) %in% rank_names[1:i])
    M <- Abar[newM.node, newM.node]
    
    ###### section for principal component analysis
    X1_PC <- data[,which(colnames(data) %in% rank_names[1:i])]
    #print(X1_PC)
    # calculate individual protein correlation
    protein_correlation <- cor(X1_PC, Pheno)
    protein_correlation_data <- data.frame(name = colnames(X1_PC), correlation = protein_correlation)
    # PCA function
    if (method == 'PCA')
    {
      # run pca
      pca_x1_summary <- prcomp(X1_PC)
      # Obtain summary
      summary_result <- summary(pca_x1_summary)
      # Extract the first three PC scores
      pca_x1_pc1 <- data.frame(pc1 = pca_x1_summary$x[,1], 
                               pc2 = pca_x1_summary$x[,2],
                               pc3 = pca_x1_summary$x[,3],
                               y = Pheno)
    }else if(method == 'NetSHy'){
      pca_x1_summary <- hybrid_score(X1_PC, M, 
                                     npc = 3, is_alpha = FALSE)
      # Extract the first three PC scores
      pca_x1_pc1 <- data.frame(pc1 = pca_x1_summary[[1]][,1], 
                               pc2 = pca_x1_summary[[1]][,2],
                               pc3 = pca_x1_summary[[1]][,3],
                               y = Pheno)
      
    }else{
      stop('Either method of PCA or NetSHy should be provided.')
    } 
    
    
    pc_correlation <- cor(pca_x1_pc1[,1:3], pca_x1_pc1[,4])
    summary_correlation <- c(summary_correlation, pc_correlation[1])
    #cor_index <- which.max(abs(cor(pca_x1_pc1, Y)))
    if (i == min_mod_size)
      score_frame <- data.frame(pca_x1_pc1$pc1)
    else
      score_frame = cbind(score_frame, pca_x1_pc1$pc1)
  }
  
  corr_pheno <- abs(cor(score_frame, Pheno))
  
  candidate_size_1 <- min(which(corr_pheno > (0.9 * max(corr_pheno))))
  cormat <-  abs(round(x = cor(score_frame[,candidate_size_1:(max_mod_size - min_mod_size + 1)]), digits = 2))
  candidate_size_2 <- max(which(cormat[,1] > 0.8))
  print(candidate_size_2)
  mod_size <- candidate_size_2 + candidate_size_1 - 2 + min_mod_size
  
  
  # finally calculate the optimal network size 
  #  mod_size <- selection(cormat[,1], summary_correlation,cor_cut = 0.8, default_size = min_mod_size, network_preference = network_preference)
  print(mod_size)
  # obtain optimal result
  newM.node <- which(colnames(Abar) %in% rank_names[1:mod_size])
  sub_type <- type[newM.node]
  M <- Abar[newM.node, newM.node]
  
  # print(corr_pheno)
  
  
  ###### section for principal component analysis
  X1_PC <- data[,which(colnames(data) %in% rank_names[1:mod_size])]
  # calculate individual protein correlation
  protein_correlation <- cor(X1_PC, Pheno)
  protein_correlation_test <- rep(0, ncol(X1_PC))
  for (i in 1:ncol(X1_PC))
  {
    protein_correlation_test[i] <- cor.test(X1_PC[,i], Pheno)$p.value
  }
  omics_correlation_data <- data.frame(name = colnames(X1_PC), correlation = protein_correlation, p = protein_correlation_test)
  if (method == 'PCA')
  {
    # run pca
    pca_x1_summary <- prcomp(X1_PC)
    # Obtain summary
    summary_result <- summary(pca_x1_summary)
    # Extract the first three PC scores
    pca_x1_pc1 <- data.frame(pc1 = pca_x1_summary$x[,1], 
                             pc2 = pca_x1_summary$x[,2],
                             pc3 = pca_x1_summary$x[,3],
                             y = Pheno)
    pc_loading <- summary_result[["rotation"]]
  }else if(method == 'NetSHy'){
    pca_x1_summary <- hybrid_score(X1_PC, M, 
                                   npc = 3, is_alpha = FALSE)
    # Extract the first three PC scores
    pca_x1_pc1 <- data.frame(pc1 = pca_x1_summary[[1]][,1], 
                             pc2 = pca_x1_summary[[1]][,2],
                             pc3 = pca_x1_summary[[1]][,3],
                             y = Pheno)
    pc_loading <- pca_x1_summary[[3]]
    
  } 
  pc_correlation <- cor(pca_x1_pc1[,1:3], pca_x1_pc1[,4])
  summary_correlation <- c(summary_correlation, pc_correlation[1])
  
  correlation_sub <- CorrMatrix[newM.node, newM.node]
  
  
  
  
  
  # Save all the cc results into a same data file
  save(M = M, pca_score = pca_x1_pc1,pc_loading, rank_value = rank_value, pc_correlation, omics_correlation_data,
       mod_size, sub_type, summary_correlation, correlation_sub,
       file = paste0(saving_dir, "/size_", mod_size,"_net_",ModuleIdx,".Rdata"))
  
  cat(paste0('The final network size is: ', nrow(M), ' with PC1 correlation w.r.t. phenotype to be: ', round(pc_correlation[1], 3)))
  #return(c(nrow(M),pc_correlation))
}

