# set up covariance tensor 
library(dplyr)
###### 4-way covariance tensor
COV_tensor <- array(0, dim = c(nrow(X1), nrow(X2), nrow(X3), nrow(Y1)))
for (i in 1:dim(X1)[2])
{
  C <- outer(X1[,i], X2[,i])
  C <- outer(C, X3[,i])
  C <- outer(C, Y1[,i])
  COV_tensor <- COV_tensor + C
  print(i)
}
COV_tensor <- abs(COV_tensor/(ncol(X1)))


###### Covariance Tensor 1 and 4
COV_14 <- array(rep(0, nrow(X1) * nrow (Y1)), dim = c(nrow(X1),  nrow(Y1)))
for (i in 1:dim(X1)[2])
{
  C <- outer(X1[,i], Y1[,i])
  COV_14 <- COV_14 + C
}
COV_14 <- abs(COV_14/(ncol(X1)))


###### Covariance Tensor 2 and 4
COV_24 <- array(rep(0, nrow(X2) * nrow (Y1)), dim = c(nrow(X2),  nrow(Y1)))
for (i in 1:dim(X2)[2])
{
  C <- outer(X2[,i], Y1[,i])
  COV_24 <- COV_24 + C
}
COV_24 <- abs(COV_24/(ncol(X2)))

###### Covariance Tensor 3 and 4
COV_34 <- array(rep(0, nrow(X3) * nrow (Y1)), dim = c(nrow(X3),  nrow(Y1)))
for (i in 1:dim(X1)[2])
{
  C <- outer(X3[,i], Y1[,i])
  COV_34 <- COV_34 + C
}
COV_34 <- abs(COV_34/(ncol(X3)))




###### Covariance Tensor 1, 2 and 4
COV_124 <- array(rep(0, nrow(X1) * nrow(X2) * nrow (Y1)), dim = c(nrow(X1),nrow(X2),nrow(Y1)))
for (i in 1:dim(X2)[2])
{
  C <- outer(X1[,i], X2[,i])
  C <- outer(C, Y1[,i])
  COV_124 <- COV_124 + C
}
COV_124 <- abs(COV_124/(ncol(Y1)))


COV_tensor <- as.tensor(COV_tensor)
COV_14 <- as.tensor(COV_14)
COV_24 <- as.tensor(COV_24)
COV_34 <- as.tensor(COV_34)
COV_124 <- as.tensor(COV_124)
##############################################################################################

# Step 2: calculate biased sampling vector

###### Set all-one vector
x1_one <- t(as.matrix(rep(1, nrow(X1))))
x2_one <- t(as.matrix(rep(1, nrow(X2))))
x3_one <- t(as.matrix(rep(1, nrow(X3))))



###### Define Probability Vector
samp_prob_1 <-  ttl(COV_tensor,list(x2_one, x3_one),ms=2:3)@data[,1,1,1]/sum(ttl(COV_tensor,list(x2_one, x3_one),ms=2:3)@data[,1,1,1]) +
  ttl(COV_124,list(x2_one),ms=2)@data[,1,1]/sum(ttl(COV_124,list(x2_one),ms=2)@data[,1,1]) + 
  COV_14@data[,1]/sum(COV_14@data[,1])


samp_prob_2 <-  ttl(COV_tensor,list(x1_one, x3_one),ms=c(1,3))@data[1,,1,1]/sum(ttl(COV_tensor,list(x1_one, x3_one),ms=c(1,3))@data[1,,1,1]) +
  ttl(COV_124,list(x1_one),ms=1)@data[1,,1]/sum(ttl(COV_124,list(x1_one),ms=1)@data[1,,1]) + 
  COV_24@data[,1]/sum(COV_24@data[,1])

samp_prob_3 <-  ttl(COV_tensor,list(x2_one, x3_one),ms=2:3)@data[1,1,,1]/sum(ttl(COV_tensor,list(x2_one, x3_one),ms=2:3)@data[1,1,,1]) +
  COV_34@data[,1]/sum(COV_34@data[,1])

samp_prob_1 <- abs(COV_14@data)^2
samp_prob_2 <- abs(COV_24@data)^2
samp_prob_3 <- abs(COV_34@data)^2

cor(t(X1), t(Y1))

###### cache covariance tensor
COV_tensor <- NULL
COV_14 <- NULL
COV_24 <- NULL
COV_34 <- NULL
COV_124 <- NULL

##############################################################################################

# step 3: GTCCA algorithm


###### Define initial parameter
common_percent <- 0.2
additional_percent <- 0.1
num_ite <- 10

######: GTCCA built-in function

turbo_gtcca <- function(num_ite = 5, common_percent = 0.1, additional_percent = 0.03, X1, X2, X3, Y1, num_modules = 3,
                        seed = 123456, prob_vec)
{
  # Define Empty canonical weight array to store vectors
  CC_X1 <- array(0, dim = c(nrow(X1), num_modules))
  CC_X2 <- array(0, dim = c(nrow(X2), num_modules))
  CC_X3 <- array(0, dim = c(nrow(X3), num_modules))
  
  # Create placeholder
  Temp_X1 <- array(0, dim = c(nrow(X1), num_modules))
  Temp_X2 <- array(0, dim = c(nrow(X2), num_modules))
  Temp_X3 <- array(0, dim = c(nrow(X3), num_modules))
  
  # random common index
  #common_index_x1 <- sample(1:nrow(X1), ceiling(nrow(X1) * common_percent),replace = FALSE, prob = prob_vec[[1]])
  #common_index_x2 <- sample(1:nrow(X2), ceiling(nrow(X2) * common_percent),replace = FALSE, prob = prob_vec[[2]])
  #common_index_x3 <- sample(1:nrow(X3), ceiling(nrow(X3) * common_percent),replace = FALSE, prob = prob_vec[[3]])
  
  # top common index
  num_common_1 = ceiling(nrow(X1) * common_percent)
  num_common_2 = ceiling(nrow(X2) * common_percent)
  num_common_3 = ceiling(nrow(X3) * common_percent)
  common_index_x1 <- order(prob_vec[[1]], decreasing = TRUE)[1:num_common_1]
  common_index_x2 <- order(prob_vec[[2]], decreasing = TRUE)[1:num_common_2]
  common_index_x3 <- order(prob_vec[[3]], decreasing = TRUE)[1:num_common_3]
  
  # randomly select additional samples
  complete_index_x1 <- 1:nrow(X1)
  complete_index_x2 <- 1:nrow(X2)
  complele_index_x3 <- 1:nrow(X3)
  
  overall_merge_index_x1 <- common_index_x1
  overall_merge_index_x2 <- common_index_x2
  overall_merge_index_x3 <- common_index_x3
  
  
  # Store subsamples 
  X1_sub <- array(0, dim = c(nrow(X1), num_modules, num_ite))
  X2_sub <- array(0, dim = c(nrow(X2), num_modules, num_ite))
  X3_sub <- array(0, dim = c(nrow(X3), num_modules, num_ite))
  
  for (ite in 1:num_ite)
  {
    set.seed(ite)
    additional_index_x1 <- sample(complete_index_x1[-common_index_x1],ceiling(nrow(X1) * additional_percent), replace = FALSE, prob = prob_vec[[1]][-common_index_x1])
    additional_index_x2 <- sample(complete_index_x2[-common_index_x2],ceiling(nrow(X2) * additional_percent), replace = FALSE, prob = prob_vec[[2]][-common_index_x2])
    additional_index_x3 <- sample(complele_index_x3[-common_index_x3],ceiling(nrow(X3) * additional_percent), replace = FALSE, prob = prob_vec[[3]][-common_index_x3])
    #if (ite == num_ite)
    #{
     # print('hahaha')
      #num_additional_1 = ceiling(nrow(X1) * additional_percent)
      #num_additional_2 = ceiling(nrow(X2) * additional_percent)
      #num_additional_3 = ceiling(nrow(X3) * additional_percent)
      #a1 <- 1:nrow(X1)
      #a2 <- 1:nrow(X2)
      #a3 <- 1:nrow(X3)
      #additional_index_x1 <- a1[-overall_merge_index_x1][order(prob_vec[[1]][-overall_merge_index_x1], decreasing = TRUE)][1:num_additional_1]
      #additional_index_x2 <- a2[-overall_merge_index_x2][order(prob_vec[[2]][-overall_merge_index_x2], decreasing = TRUE)][1:num_additional_2]
      #additional_index_x3 <- a3[-overall_merge_index_x3][order(prob_vec[[3]][-overall_merge_index_x3], decreasing = TRUE)][1:num_additional_3]
    #}
    # merge index 
    merge_index_x1 <- unique(sort(c(common_index_x1, additional_index_x1)))
    merge_index_x2 <- unique(sort(c(common_index_x2, additional_index_x2)))
    merge_index_x3 <- unique(sort(c(common_index_x3, additional_index_x3)))
    
    # merge index 
    overall_merge_index_x1 <- unique(sort(c(overall_merge_index_x1, additional_index_x1)))
    overall_merge_index_x2 <- unique(sort(c(overall_merge_index_x2, additional_index_x2)))
    overall_merge_index_x3 <- unique(sort(c(overall_merge_index_x3, additional_index_x3)))
    print(overall_merge_index_x1)
    # subset data
    X1_subset <- X1[merge_index_x1, ]
    X2_subset <- X2[merge_index_x2, ]
    X3_subset <- X3[merge_index_x3, ]
    
    
    
    tryCatch(
      # starting with the first iteration
      if (ite == 1)
      {
        # extract canonical weight
        result_weight <- gtcca(X1_subset, X2_subset, X3_subset, Y1, num_modules)
        
        print(dim(result_weight[[1]]))
        print(dim(CC_X1[merge_index_x1,]))
        print(result_weight[[1]])
        
        # distribute canonical weight
        CC_X1[merge_index_x1,] <- result_weight[[1]]
        CC_X2[merge_index_x2,] <- result_weight[[2]]
        CC_X3[merge_index_x3,] <- result_weight[[3]]
        
        
        
        
        # normalization by common index
        CC_X1 <- apply(CC_X1, 2, function(x){x/pracma::Norm(x[common_index_x1])})
        CC_X2 <- apply(CC_X2, 2, function(x){x/pracma::Norm(x[common_index_x2])})
        CC_X3 <- apply(CC_X3, 2, function(x){x/pracma::Norm(x[common_index_x3])})
        
        X1_sub[,,1] <- CC_X1
        X2_sub[,,1] <- CC_X2
        X3_sub[,,1] <- CC_X3
        
        # store the result for the first iteration
        initialize_x1 <- CC_X1
        initialize_x2 <- CC_X2
        initialize_x3 <- CC_X3
        
        # store initial sub-sample index for comparison in later iterations
        initial_merge_index_x1 <- merge_index_x1
        initial_merge_index_x2 <- merge_index_x2
        initial_merge_index_x3 <- merge_index_x3
        
        
        
        
        
      }
      else{
        # extract canonical weight
        result_weight <- gtcca(X1_subset, X2_subset, X3_subset, Y1, num_modules)
        
        # distribute canonical weight
        Temp_X1[merge_index_x1,] <- result_weight[[1]]
        Temp_X2[merge_index_x2,] <- result_weight[[2]]
        Temp_X3[merge_index_x3,] <- result_weight[[3]]
        
        # normalization by common index
        Temp_X1 <- apply(Temp_X1, 2, function(x){x/pracma::Norm(x[common_index_x1])})
        Temp_X2 <- apply(Temp_X2, 2, function(x){x/pracma::Norm(x[common_index_x2])})
        Temp_X3 <- apply(Temp_X3, 2, function(x){x/pracma::Norm(x[common_index_x3])})
        
        # keep track of which correspondence hasn't been assigned yet
        track_index <- 1:num_modules
        # finding solution correspondence by the first view
        for (solution in 1:num_modules)
        { 
          # finding the max correlation
          max_cor_index <- which.max(abs(cor(Temp_X1[common_index_x1,solution],initialize_x1[common_index_x1,track_index])))
          print(max(abs(cor(Temp_X1[common_index_x1,solution],initialize_x1[common_index_x1,track_index]))))
          print(max_cor_index)
          # stitch algorithm
          if(cor(Temp_X1[common_index_x1,solution],initialize_x1[common_index_x1,max_cor_index]) < 0)
            Temp_X1[,solution] <- -Temp_X1[,solution]
          
          # stitch unseen weight
          CC_X1[merge_index_x1 %!in% initial_merge_index_x1, track_index[max_cor_index]] <- Temp_X1[merge_index_x1 %!in% initial_merge_index_x1, solution]
          CC_X2[merge_index_x2 %!in% initial_merge_index_x2, track_index[max_cor_index]] <- Temp_X2[merge_index_x2 %!in% initial_merge_index_x2, solution]
          CC_X3[merge_index_x3 %!in% initial_merge_index_x3, track_index[max_cor_index]] <- Temp_X3[merge_index_x3 %!in% initial_merge_index_x3, solution]
          print(ite)
          X1_sub[merge_index_x1,track_index[max_cor_index],ite] <- Temp_X1[merge_index_x1, solution]
          X2_sub[merge_index_x2,track_index[max_cor_index],ite] <- Temp_X2[merge_index_x2, solution]
          X3_sub[merge_index_x3,track_index[max_cor_index],ite] <- Temp_X3[merge_index_x3, solution]
          
          
          
          # eliminate one correspondence  
          track_index <- track_index[!track_index %in% track_index[max_cor_index]]  
          
          
          
          
        }
        print('new indexes based on new iteration are: ')
        print(merge_index_x1[merge_index_x1 %!in% initial_merge_index_x1])
        
        # update initial merge index
        initial_merge_index_x1 <- union(initial_merge_index_x1, merge_index_x1)
        initial_merge_index_x2 <- union(initial_merge_index_x2, merge_index_x2)
        initial_merge_index_x3 <- union(initial_merge_index_x3, merge_index_x3)
        
        
      })
    
    print(paste0('The ', ite, 'th iteration is done!'))
  }  
  return(list(CC_X1, CC_X2, CC_X3, X1_sub, X2_sub, X3_sub))
}





###################################################### GTCCA Source Function


gtcca <- function(X1, X2, X3, Y1, num_modules)
{
  ###### 4-way covariance tensor
  COV_tensor <- array(0, dim = c(nrow(X1), nrow(X2), nrow(X3), nrow(Y1)))
  for (i in 1:dim(X1)[2])
  {
    C <- outer(X1[,i], X2[,i])
    C <- outer(C, X3[,i])
    C <- outer(C, Y1[,i])
    COV_tensor <- COV_tensor + C
  }
  COV_tensor <- abs(COV_tensor/(ncol(X1)))
  
  
  ###### Covariance Tensor 1 and 4
  COV_14 <- array(rep(0, nrow(X1) * nrow (Y1)), dim = c(nrow(X1),  nrow(Y1)))
  for (i in 1:dim(X1)[2])
  {
    C <- outer(X1[,i], Y1[,i])
    COV_14 <- COV_14 + C
  }
  COV_14 <- abs(COV_14/(ncol(X1)))
  
  
  ###### Covariance Tensor 2 and 4
  COV_24 <- array(rep(0, nrow(X2) * nrow (Y1)), dim = c(nrow(X2),  nrow(Y1)))
  for (i in 1:dim(X2)[2])
  {
    C <- outer(X2[,i], Y1[,i])
    COV_24 <- COV_24 + C
  }
  COV_24 <- abs(COV_24/(ncol(X2)))
  
  ###### Covariance Tensor 3 and 4
  COV_34 <- array(rep(0, nrow(X3) * nrow (Y1)), dim = c(nrow(X3),  nrow(Y1)))
  for (i in 1:dim(X1)[2])
  {
    C <- outer(X3[,i], Y1[,i])
    COV_34 <- COV_34 + C
  }
  COV_34 <- abs(COV_34/(ncol(X3)))
  
  
  
  
  ###### Covariance Tensor 1, 2 and 4
  COV_124 <- array(rep(0, nrow(X1) * nrow(X2) * nrow (Y1)), dim = c(nrow(X1),nrow(X2),nrow(Y1)))
  for (i in 1:dim(X2)[2])
  {
    C <- outer(X1[,i], X2[,i])
    C <- outer(C, Y1[,i])
    COV_124 <- COV_124 + C
  }
  COV_124 <- abs(COV_124/(ncol(Y1)))
  
  # convert into tensor format
  COV_tensor <- as.tensor(COV_tensor)
  COV_14 <- as.tensor(COV_14)
  COV_24 <- as.tensor(COV_24)
  COV_34 <- as.tensor(COV_34)
  COV_124 <- as.tensor(COV_124)
  
  
  
  
  cc_tensor <- rep(NA, num_modules)
  cc_14 <- rep(NA, num_modules)
  cc_24 <- rep(NA, num_modules)
  cc_34 <- rep(NA, num_modules)
  cc_124 <- rep(NA, num_modules)
  cc_ss <- rep(NA, num_modules)
  CCcoef <- c(1,1,1,1,1)
  
  
  
  assign('gene_vector',matrix(nrow = nrow(X1)))
  assign('protein_vector',matrix(nrow = nrow(X2)))
  assign('metabolite_vector',matrix(nrow = nrow(X3)))
  
  reconstruct_tensor <- array(0, dim = c(nrow(X1), nrow(X2), nrow(X3), nrow(Y1)))
  reconstruct_14 <- array(0, dim = c(nrow(X1), nrow(Y1)))
  reconstruct_24 <- array(0, dim = c(nrow(X2), nrow(Y1)))
  reconstruct_34 <- array(0, dim = c(nrow(X3), nrow(Y1)))
  reconstruct_124 <- array(0, dim = c(nrow(X1), nrow(X2), nrow(Y1)))
  
  
  cov_tensor <- COV_tensor
  cov_14 <- COV_14
  cov_24 <- COV_24
  cov_34 <- COV_34
  cov_124 <- COV_124
  
  
  
  
  cc_tensor <- rep(NA, num_modules)
  cc_14 <- rep(NA, num_modules)
  cc_24 <- rep(NA, num_modules)
  cc_34 <- rep(NA, num_modules)
  cc_124 <- rep(NA, num_modules)
  cc_ss <- rep(NA, num_modules)
  CCcoef <- c(1,1,1,1,1)
  
  
  
  assign('metabolite_vector',matrix(nrow = nrow(X1)))
  assign('protein_vector',matrix(nrow = nrow(X2)))
  assign('gene_vector',matrix(nrow = nrow(X3)))
  assign('pheno_vector',matrix(nrow = nrow(Y1)))
  
  reconstruct_tensor <- array(0, dim = c(nrow(X1), nrow(X2), nrow(X3), nrow(Y1)))
  reconstruct_14 <- array(0, dim = c(nrow(X1), nrow(Y1)))
  reconstruct_24 <- array(0, dim = c(nrow(X2), nrow(Y1)))
  reconstruct_34 <- array(0, dim = c(nrow(X3), nrow(Y1)))
  reconstruct_124 <- array(0, dim = c(nrow(X1), nrow(X2), nrow(Y1)))
  
  
  for (j in 1:num_modules)
  {  
    
    r <- 1
    
    u_1 <- matrix(rnorm(nrow(X1) * r, mean = 0 , sd = 1), nrow = nrow(X1))
    u_1 <- apply(u_1, 2, function(x){x/pracma::Norm(x)})
    u_2 <- matrix(rnorm(nrow(X2) * r, mean = 0 , sd = 1), nrow = nrow(X2))
    u_2 <- apply(u_2, 2, function(x){x/pracma::Norm(x)})
    u_3 <- matrix(rnorm(nrow(X3) * r, mean = 0 , sd = 1), nrow = nrow(X3))
    u_3 <- apply(u_3, 2, function(x){x/pracma::Norm(x)})
    u_4 <- matrix(rnorm(nrow(Y1) * r, mean = 0 , sd = 1), nrow = nrow(Y1))
    u_4 <- apply(u_4, 2, function(x){x/pracma::Norm(x)})
    
    
    
    vec <- c(as.vector(u_1),as.vector(u_2),as.vector(u_3),as.vector(u_4))
    vec <- as.numeric(vec)
    
    superdiagonal_tensor <- function(num_modes,len,elements=1L){
      modes <- rep(len,num_modes)
      arr <- array(0, dim = modes)
      if(length(elements)==1) elements <- rep(elements,len)
      for (i in 1:len){
        txt <- paste("arr[",paste(rep("i", num_modes),collapse=","),"] <- ", elements[i],sep="")
        eval(parse(text=txt))
      }
      as.tensor(arr)
    }
    
    index_cumsum <- cumsum(c(nrow(X1), nrow(X2), nrow(X3), nrow(Y1)))
    
    opt_fun <- function(vec)
    {
      u_1 <- matrix(vec[1:(index_cumsum[1]*r)], ncol = r)
      u_2 <- matrix(vec[(index_cumsum[1]*r+1):(index_cumsum[2]*r)], ncol = r)
      u_3 <- matrix(vec[(index_cumsum[2] * r + 1) :(index_cumsum[3] * r)], ncol = r)
      u_4 <- matrix(vec[(index_cumsum[3] * r + 1):(index_cumsum[4] * r)], ncol = r)
      
      sup_lambda <- superdiagonal_tensor(4, r, 1)
      recon <-  ttl(sup_lambda,list(u_1, u_2, u_3, u_4),ms=1:4)
      
      sup_rho <- superdiagonal_tensor(2, r, 1)
      recon_14 <- ttl(sup_rho,list(u_1, u_4),ms=1:2)
      
      sup_rho1 <- superdiagonal_tensor(2, r, 1)
      recon_24 <- ttl(sup_rho1,list(u_2, u_4),ms=1:2)
      
      sup_rho2 <- superdiagonal_tensor(2, r, 1)
      recon_34 <- ttl(sup_rho2,list(u_3, u_4),ms=1:2)
      
      sup_rho3 <- superdiagonal_tensor(3, r, 1)
      recon_124 <- ttl(sup_rho3,list(u_1,u_2, u_4),ms=1:3)
      
      
      norm_u1 <- (apply(u_1, 2, function(x){pracma::Norm(x)})-1)^2
      norm_u2 <- (apply(u_2, 2, function(x){pracma::Norm(x)})-1)^2
      norm_u3 <- (apply(u_3, 2, function(x){pracma::Norm(x)})-1)^2
      norm_u4 <- (apply(u_4, 2, function(x){pracma::Norm(x)})-1)^2
      
      #recon <- as.tensor(lambda * outer(outer(u_1, u_2), u_3))
      #recon_14 <- as.tensor(rho * outer(u_1, u_4))
      #recon_34 <- as.tensor(rho_2 * outer(u_3, u_4))
      #recon_24 <- as.tensor(rho_1 * outer(u_2, u_4))
      return(CCcoef[1] * fnorm(COV_tensor - recon)^2 + CCcoef[2] * fnorm(COV_14 - recon_14)^2 + CCcoef[3] * fnorm(COV_24 - recon_24)^2 +  CCcoef[4] * fnorm(COV_34 - recon_34)^2 + CCcoef[5] * fnorm(COV_124 - recon_124)^2)
    }
    
    der_fun <- function(vec)
    {
      u_1 <- matrix(vec[1:(index_cumsum[1]*r)], ncol = r)
      u_2 <- matrix(vec[(index_cumsum[1]*r+1):(index_cumsum[2]*r)], ncol = r)
      u_3 <- matrix(vec[(index_cumsum[2] * r + 1) :(index_cumsum[3] * r)], ncol = r)
      u_4 <- matrix(vec[(index_cumsum[3] * r + 1):(index_cumsum[4] * r)], ncol = r)
      
      sup_lambda <- superdiagonal_tensor(4, r, 1)
      recon <-  ttl(sup_lambda,list(u_1, u_2, u_3, u_4),ms=1:4)
      
      sup_rho <- superdiagonal_tensor(2, r, 1)
      recon_14 <- ttl(sup_rho,list(u_1, u_4),ms=1:2)
      
      sup_rho1 <- superdiagonal_tensor(2, r, 1)
      recon_24 <- ttl(sup_rho1,list(u_2, u_4),ms=1:2)
      
      sup_rho2 <- superdiagonal_tensor(2, r, 1)
      recon_34 <- ttl(sup_rho2,list(u_3, u_4),ms=1:2)
      
      sup_rho3 <- superdiagonal_tensor(3, r, 1)
      recon_124 <- ttl(sup_rho3,list(u_1, u_2, u_4),ms=1:3)
      
      
      norm_1 <- u_1 - apply(u_1, 2, function(x){x/pracma::Norm(x)})
      norm_2 <- u_2 - apply(u_2, 2, function(x){x/pracma::Norm(x)})
      norm_3 <- u_3 - apply(u_3, 2, function(x){x/pracma::Norm(x)})
      norm_4 <- u_4 - apply(u_4, 2, function(x){x/pracma::Norm(x)})
      
      
      
      
      vec <- c(
        as.vector(CCcoef[1] * (k_unfold(recon, 1)@data - k_unfold(COV_tensor, 1)@data) %*% (khatri_rao_list(list(as.matrix(u_4), as.matrix(u_3),as.matrix(u_2)))) + CCcoef[2] * (u_1 %*% t(u_4) - COV_14@data) %*% u_4 + CCcoef[5] * (k_unfold(recon_124, 1)@data - k_unfold(COV_124, 1)@data) %*% (khatri_rao_list(list(as.matrix(u_4), as.matrix(u_2))))),
        as.vector(CCcoef[1] * (k_unfold(recon, 2)@data - k_unfold(COV_tensor, 2)@data) %*% (khatri_rao_list(list(as.matrix(u_4), as.matrix(u_3),as.matrix(u_1)))) + CCcoef[3] * (u_2 %*% t(u_4) - COV_24@data) %*% u_4 + CCcoef[5] * (k_unfold(recon_124, 2)@data - k_unfold(COV_124, 2)@data) %*% (khatri_rao_list(list(as.matrix(u_4), as.matrix(u_1))))),
        as.vector(CCcoef[1] * (k_unfold(recon, 3)@data - k_unfold(COV_tensor, 3)@data) %*% (khatri_rao_list(list(as.matrix(u_4), as.matrix(u_2),as.matrix(u_1)))) + CCcoef[4] * (u_3 %*% t(u_4) - COV_34@data) %*% u_4),
        as.vector(CCcoef[1] * (k_unfold(recon, 4)@data - k_unfold(COV_tensor, 4)@data) %*% (khatri_rao_list(list(as.matrix(u_3), as.matrix(u_2),as.matrix(u_1)))) + CCcoef[5] * (k_unfold(recon_124, 3)@data - k_unfold(COV_124, 3)@data) %*% (khatri_rao_list(list(as.matrix(u_2), as.matrix(u_1))))+ CCcoef[2] * t(u_1 %*% t(u_4) - COV_14@data) %*% u_1 + CCcoef[3] * t(u_2 %*% t(u_4) - COV_24@data) %*% u_2 + CCcoef[4] * t(u_3 %*% t(u_4) - COV_34@data) %*% u_3)
      )
      return(as.numeric(vec))
    }
    
    
    system.time(
      result <- Rcgmin::Rcgmin(par = vec, fn = opt_fun, gr = der_fun,  control = list(eps = 1e-07, maxit = 150))
    )
    
    
    u_1 <- matrix(result$par[1:(index_cumsum[1]*r)], ncol = r)
    u_2 <- matrix(result$par[(index_cumsum[1]*r+1):(index_cumsum[2]*r)], ncol = r)
    u_3 <- matrix(result$par[(index_cumsum[2] * r + 1) :(index_cumsum[3] * r)], ncol = r)
    u_4 <- matrix(result$par[(index_cumsum[3] * r + 1):(index_cumsum[4] * r)], ncol = r)
    
    
    
    
    gene_vector <- cbind(gene_vector,u_1)
    protein_vector <- cbind(protein_vector,u_2)
    metabolite_vector <- cbind(metabolite_vector,u_3)
    pheno_vector <- cbind(pheno_vector,u_4)
    
    
    
    
    sup_lambda <- superdiagonal_tensor(4, r, 1)
    recon <-  ttl(sup_lambda,list(u_1, u_2, u_3, u_4),ms=1:4)
    sup_rho <- superdiagonal_tensor(2, r, 1)
    recon_14 <- ttl(sup_rho,list(u_1, u_4),ms=1:2)
    sup_rho1 <- superdiagonal_tensor(2, r, 1)
    recon_24 <- ttl(sup_rho1,list(u_2, u_4),ms=1:2)
    sup_rho2 <- superdiagonal_tensor(2, r, 1)
    recon_34 <- ttl(sup_rho2,list(u_3, u_4),ms=1:2)
    sup_rho3 <- superdiagonal_tensor(3, r, 1)
    recon_124 <- ttl(sup_rho3,list(u_1,u_2, u_4),ms=1:3)
    
    
    reconstruct_tensor <- reconstruct_tensor + recon@data
    reconstruct_14 <- reconstruct_14 + recon_14@data
    reconstruct_24 <- reconstruct_24 + recon_24@data
    reconstruct_34 <- reconstruct_34 + recon_34@data
    reconstruct_124 <- reconstruct_124 + recon_124@data
    
    
    
    
    
    
    COV_tensor <- COV_tensor - recon
    COV_14 <- COV_14 - recon_14
    COV_24 <- COV_24 - recon_24
    COV_34 <- COV_34 - recon_34
    COV_124 <- COV_124 - recon_124
    
    
    
    
    cov_tensor_var <- 1- (fnorm(COV_tensor)/fnorm(cov_tensor))
    cov_14_var <- 1- (fnorm(COV_14)/fnorm(cov_14))
    cov_24_var <- 1- (fnorm(COV_24)/fnorm(cov_24))
    cov_34_var <- 1- (fnorm(COV_34)/fnorm(cov_34))
    cov_124_var <- 1- (fnorm(COV_124)/fnorm(cov_124))
    
    
  } 
  return(list(gene_vector[,-1], protein_vector[,-1], metabolite_vector[,-1]))
  
}

# define not in for use
'%!in%' <- function(x,y)!('%in%'(x,y))


# start the algorithm
result <- turbo_gtcca(X1 = X1, X2 = X2, X3 = X3, Y1 = Y1, prob_vec = list(samp_prob_1, samp_prob_2, samp_prob_3))
result[[1]]
result[[1]] <- apply(result[[1]],2, function(x){return(x/pracma::Norm(x))})
agg_gtcca <- apply(rbind(abs(result[[1]]),abs(result[[2]]),abs(result[[3]])),1, function(x) sum(x))


a <- 1:1000
dim(result[[4]])
plot(abs(result[[1]][,1])~a,type = 'h')
empty <- matrix(0, nrow = 3000, ncol = 5)
for (i in 1:5)
{
  empty[,i] <- apply(rbind(abs(result[[4]][,,i]),abs(result[[5]][,,i]),abs(result[[6]][,,i])),1, function(x) sum(x))
}

agg_gtcca <- apply(empty, 1, mean)

plot(agg_result~a,type = 'h')
my_cor <- cor(t(X1), t(Y1))
plot(samp_prob_1~a,type = 'h')
plot(agg_gtcca, type = 'h')
plot(abs(result[[1]][,1]), type = 'h')
plot(abar1_nopen_weighted[,1], type = 'h')



library(mltest)
membership_x1 <- rep(0,150)
membership_x2 <- rep(0,150)
membership_x3 <- rep(0,150)
for (i in 1:150)
{
  membership_x1[i] <- which.max(abs(result[[1]][i,]))
  membership_x2[i] <- which.max(abs(result[[2]][i,]))
  membership_x3[i] <- which.max(abs(result[[3]][i,1:2]))
}
membership_x3 <- membership_x3[c(1:50,101:150)]
true_x1 <- c(rep(1,50), rep(2, 50), rep(3,50))
true_x3 <- c(rep(1,50), rep(2,50))
ml_test(membership_x1, true_x1)$accuracy
ml_test(membership_x2, true_x1)
ml_test(membership_x3, true_x3)
###
dim(result[[4]])
omics_1 <- matrix(0, nrow = 1000, ncol = 3)
for (i in 1:10)
{
  temp <- result[[4]][i,,]
  print(temp)
}


######################## TCCA
cp(tnsr, num_components = NULL, max_iter = 25, tol = 1e-05)



####################### unweighted smcca no penalty
L1 <- max(1, sqrt(nrow(X1)))
L2 <- max(1, sqrt(nrow(X2)))
L3 <- max(1, sqrt(nrow(X3)))
L4 <- max(1, sqrt(nrow(Y1)))

Ws <- myMultiCCA(xlist = list(t(X1),t(X2),t(X3),t(Y1)), penalty = c(L1,L2,L3,L4), type = 'standard',
                 ncomponents = 3, trace = TRUE,standardize = TRUE, CCcoef = c(1,1,1,1,1,1))
abar1_nopen_unweighted <- as.matrix(abs(Ws$ws[[1]]))
abar2_nopen_unweighted <- as.matrix(abs(Ws$ws[[2]]))
abar3_nopen_unweighted <- as.matrix(abs(Ws$ws[[3]]))
agg_smcca_nopen_unweighted <- apply(rbind(abar1_nopen_unweighted, abar2_nopen_unweighted, abar3_nopen_unweighted),1, function(x) sum(x))
plot(agg_smcca_nopen_unweighted, type = 'h')
membership_x1 <- rep(0,150)
membership_x2 <- rep(0,150)
membership_x3 <- rep(0,150)
for (i in 1:150)
{
  membership_x1[i] <- which.max(abar1_nopen_unweighted[i,])
  membership_x2[i] <- which.max(abar2_nopen_unweighted[i,])
  membership_x3[i] <- which.max(abar3_nopen_unweighted[i,1:2])
}
membership_x3 <- membership_x3[c(1:50,101:150)]
true_x1 <- c(rep(1,50), rep(3, 50), rep(2,50))
true_x3 <- c(rep(1,50), rep(2,50))
ml_test(factor(membership_x1,levels = c(1,2,3)), true_x1)
ml_test(factor(membership_x2,levels = c(1,2,3)), true_x1)
ml_test(membership_x3, true_x3)

(0.89 * 3 + 0.67 * 3 + 0.99 * 2)/8




####################### unweighted smcca penalty


L1 <- max(1, sqrt(nrow(X1))*0.2)
L2 <- max(1, sqrt(nrow(X2))*0.2)
L3 <- max(1, sqrt(nrow(X3))*0.2)
L4 <- max(1, sqrt(nrow(Y1))*0.2)


Ws <- myMultiCCA(xlist = list(t(X1),t(X2),t(X3),t(Y1)), penalty = c(L1,L2,L3,L4), type = 'standard',
                 ncomponents = 3, trace = TRUE,standardize = TRUE, CCcoef = c(1,1,1,1,1,1))
abar1_pen_unweighted <- as.matrix(abs(Ws$ws[[1]]))
abar2_pen_unweighted <- as.matrix(abs(Ws$ws[[2]]))
abar3_pen_unweighted <- as.matrix(abs(Ws$ws[[3]]))
agg_smcca_pen_unweighted <- apply(rbind(abar1_pen_unweighted, abar2_pen_unweighted, abar3_pen_unweighted),1, function(x) sum(x))
plot(agg_smcca_pen_unweighted, type = 'h')
membership_x1 <- rep(0,150)
membership_x2 <- rep(0,150)
membership_x3 <- rep(0,150)
for (i in 1:150)
{
  membership_x1[i] <- which.max(abar1_pen_unweighted[i,])
  membership_x2[i] <- which.max(abar2_pen_unweighted[i,])
  membership_x3[i] <- which.max(abar3_pen_unweighted[i,1:2])
}
membership_x3 <- membership_x3[c(1:50,101:150)]
true_x1 <- c(rep(1,50), rep(3, 50), rep(2,50))
true_x3 <- c(rep(1,50), rep(2,50))
ml_test(factor(membership_x1,levels = c(1,2,3)), true_x1)$accuracy
ml_test(factor(membership_x2,levels = c(1,2,3)), true_x1)$accuracy
ml_test(membership_x3, true_x3)$accuracy

(0.33 * 3 + 0.33 * 3 + 0.67 * 2)/8


####################### weighted smcca no penalty
L1 <- max(1, sqrt(nrow(X1)))
L2 <- max(1, sqrt(nrow(X2)))
L3 <- max(1, sqrt(nrow(X3)))
L4 <- max(1, sqrt(nrow(Y1)))

Ws <- myMultiCCA(xlist = list(t(X1),t(X2),t(X3),t(Y1)), penalty = c(L1,L2,L3,L4), type = 'standard',
                 ncomponents = 3, trace = TRUE,standardize = TRUE, CCcoef = c(1,1,5,1,5,5))
abar1_nopen_weighted <- as.matrix(abs(Ws$ws[[1]]))
abar2_nopen_weighted <- as.matrix(abs(Ws$ws[[2]]))
abar3_nopen_weighted <- as.matrix(abs(Ws$ws[[3]]))
agg_smcca_nopen_weighted <- apply(rbind(abar1_nopen_weighted, abar2_nopen_weighted, abar3_nopen_weighted),1, function(x) sum(x))
plot(agg_smcca_nopen_weighted, type = 'h')
membership_x1 <- rep(0,150)
membership_x2 <- rep(0,150)
membership_x3 <- rep(0,150)
for (i in 1:150)
{
  membership_x1[i] <- which.max(abar1_nopen_weighted[i,])
  membership_x2[i] <- which.max(abar2_nopen_weighted[i,])
  membership_x3[i] <- which.max(abar3_nopen_weighted[i,1:2])
}
membership_x3 <- membership_x3[c(1:50,101:150)]
true_x1 <- c(rep(1,50), rep(3, 50), rep(2,50))
true_x3 <- c(rep(1,50), rep(2,50))
ml_test(factor(membership_x1,levels = c(1,2,3)), true_x1)$accuracy
ml_test(factor(membership_x2,levels = c(1,2,3)), true_x1)$accuracy
ml_test(membership_x3, true_x3)$accuracy
ml_test(membership_x3, true_x3)

(0.7 * 3 + 0.67 * 3 + 0.73 * 2)/8















################################# ROC Curve
library(caret)
truth <- rep(0,3000)
truth[1:150] <- 1
truth[1001:1150] <- 1
truth[2001:2050] <- 1
truth[2101:2150] <- 1



threshold <- seq(from =  -0.01, to = 0.4, by = 0.001)
tcca_weight <- rbind(tcca_x1, tcca_x2, tcca_x3)[,1]
gtcca_weight <- rbind(abs(result[[1]]),abs(result[[2]]),abs(result[[3]]))[,1]
nopen_unweighted_weight <- rbind(abar1_nopen_unweighted, abar2_nopen_unweighted, abar3_nopen_unweighted)[,1]
pen_unweighted_weight <- rbind(abar1_pen_unweighted, abar2_pen_unweighted, abar3_pen_unweighted)[,1]
nopen_weighted_weight <- rbind(abar1_nopen_weighted, abar2_nopen_weighted, abar3_nopen_weighted)[,1]
pen_weighted_weight <- rbind(abar1_pen_weighted, abar2_pen_weighted, abar3_pen_weighted)[,1]

plot(nopen_weighted_weight, type ='h')
plot(agg_smcca_nopen_weighted, type ='h')
plot(gtcca_weight, type ='h')
plot(agg_gtcca, type = 'h')


my_roc <- function(vec, threshold)
{
  frame <- matrix(0, nrow = length(threshold), ncol = 2)
  print(frame)
  colnames(frame) <- c("sensitivity", "1-specificity")
  for (i in 1:length(threshold))
  {
    pred <- ifelse(vec > threshold[i],1,0)
    sensitivity <- sum(pred[1:150])/150
    specificity <- sum(pred[851:1000])/850
    frame[i,] <- c(sensitivity, specificity)
  }
  return(frame)
}

roc_result <- my_roc(agg_gtcca[1:1000], threshold)
plot(roc_result[,1]~roc_result[,2], type = 'l', xlim = c(0,1), ylim = c(0,1))


library(ROCR)
## Loading required package: gplots
## 
## Attaching package: 'gplots'
## The following object is masked from 'package:stats':
## 
##     lowess
# plot a ROC curve for a single prediction run
pred <- prediction(agg_gtcca[1:1000], truth[1:1000])
perf <- performance(pred,"tpr","fpr")
plot(perf,colorize=TRUE)
plot(agg_gtcca[1:1000], type = 'h')

curve1 <- roc(truth, agg_smcca_nopen_unweighted)
curve2 <- roc(truth, agg_smcca_pen_unweighted)
curve3 <- roc(truth, agg_smcca_nopen_weighted)
curve4 <- roc(truth, agg_smcca_pen_weighted)
curve5 <- roc(truth, tcca_aggregate)
curve6 <- roc(truth, agg_gtcca)
plot(curve1, print.auc = TRUE)
plot(curve2, add=TRUE, col='red', print.auc = TRUE, print.auc.y=40)
plot(curve3, add=TRUE, col='blue')
plot(curve4, add=TRUE, col='green')
plot(curve5, add=TRUE, col='green')
plot(curve6, add=TRUE, col='green')




# ROC plot 1
roc(truth, agg_smcca_nopen_unweighted, plot=TRUE, legacy.axes=FALSE, percent=TRUE, col="salmon", lwd=2, print.auc=TRUE)
plot.roc(truth, agg_smcca_pen_unweighted, percent=TRUE, col="goldenrod", lwd=2, print.auc=TRUE, add=TRUE, print.auc.y=45)
plot.roc(truth, agg_smcca_nopen_weighted, percent=TRUE, col="lightsteelblue", lwd=2, print.auc=TRUE, add=TRUE, print.auc.y=40)
plot.roc(truth, agg_smcca_pen_unweighted, percent=TRUE, col="green", lwd=2, print.auc=TRUE, add=TRUE, print.auc.y=35)
plot.roc(truth, tcca_aggregate, percent=TRUE, col="orange", lwd=2, print.auc=TRUE, add=TRUE, print.auc.y=30)
plot.roc(truth, agg_gtcca, percent=TRUE, col="red", lwd=2, print.auc=TRUE, add=TRUE, print.auc.y=25)
title(main = "ROC Curve for Signal Detection", line = 2.5)
legend(x = "right",inset = 0,
       legend = c("MCCA", "SmCCA", "Weighted MCCA", 'Weighted SmCCA','TCCA', 'GTCCA'), 
       col=c('salmon','goldenrod', 'lightsteelblue', 'green','orange', 'red'), lwd=7, cex=.7, horiz = FALSE)
example_curve <- roc(truth, nopen_weighted_weight)
auc(example_curve)

auc(curve1)
auc(curve2)
auc(curve3)
auc(curve4)
auc(curve5)
auc(curve6)

plot(curve)

ggplot(test_data, aes(date)) + 
  geom_line(aes(y = var0, colour = "var0")) + 
  geom_line(aes(y = var1, colour = "var1"))


for (i in 1:length(threshold))
{
  
}

####################### unweighted smcca penalty
L1 <- max(1, sqrt(nrow(X1))*0.2)
L2 <- max(1, sqrt(nrow(X2))*0.2)
L3 <- max(1, sqrt(nrow(X3))*0.2)
L4 <- max(1, sqrt(nrow(Y1))*0.2)

Ws <- myMultiCCA(xlist = list(t(X1),t(X2),t(X3),t(Y1)), penalty = c(L1,L2,L3,L4), type = 'standard',
                 ncomponents = 3, trace = TRUE,standardize = TRUE, CCcoef = c(1,1,5,1,5,5))
abar1_pen_weighted <- as.matrix(abs(Ws$ws[[1]]))
abar2_pen_weighted <- as.matrix(abs(Ws$ws[[2]]))
abar3_pen_weighted <- as.matrix(abs(Ws$ws[[3]]))
agg_smcca_pen_weighted <- apply(rbind(abar1_pen_weighted, abar2_pen_weighted, abar3_pen_weighted),1, function(x) sum(x))
plot(agg_smcca_pen_weighted, type = 'h')
membership_x1 <- rep(0,150)
membership_x2 <- rep(0,150)
membership_x3 <- rep(0,150)
for (i in 1:150)
{
  membership_x1[i] <- which.max(abar1_pen_weighted[i,])
  membership_x2[i] <- which.max(abar2_pen_weighted[i,])
  membership_x3[i] <- which.max(abar3_pen_weighted[i,1:2])
}
membership_x3 <- membership_x3[c(1:50,101:150)]
true_x1 <- c(rep(2,50), rep(1, 50), rep(3,50))
true_x3 <- c(rep(2,50), rep(1,50))
ml_test(factor(membership_x1,levels = c(1,2,3)), true_x1)$accuracy
ml_test(factor(membership_x2,levels = c(1,2,3)), true_x1)$accuracy
ml_test(membership_x3, true_x3)$accuracy

(0.67 * 3 + 0.65 * 3 + 0.98 * 2)/8


############################ TCCA 
load('TCCA_Weight.Rdata')
tcca_x1 <- apply(abs(check[[1]]),2, function(x){x/pracma::Norm(x)})
tcca_x2 <- apply(abs(check[[2]]),2, function(x){x/pracma::Norm(x)})
tcca_x3 <- apply(abs(check[[3]]),2, function(x){x/pracma::Norm(x)})

tcca_aggregate <- apply(rbind(tcca_x1, tcca_x2, tcca_x3),1, function(x) sum(x))
plot(tcca_aggregate, type = 'h')

membership_x1 <- rep(0,150)
membership_x2 <- rep(0,150)
membership_x3 <- rep(0,150)
for (i in 1:150)
{
  membership_x1[i] <- which.max(tcca_x1[i,])
  membership_x2[i] <- which.max(tcca_x2[i,])
  membership_x3[i] <- which.max(tcca_x3[i,1:2])
}
membership_x3 <- membership_x3[c(1:50,101:150)]
true_x1 <- c(rep(1,50), rep(3, 50), rep(2,50))
true_x3 <- c(rep(1,50), rep(2,50))
ml_test(factor(membership_x1,levels = c(1,2,3)), true_x1)$accuracy
ml_test(factor(membership_x2,levels = c(1,2,3)), true_x1)$accuracy
ml_test(membership_x3, true_x3)$accuracy

(3 * 0.35 + 3 * 0.44 + 2 * 0.78)/8



####################### one-omics
X_big <- rbind(X1,X2,X3)
L1 <- max(1, sqrt(nrow(X_big)) * 0.2)
L2 <- max(1, sqrt(nrow(Y1)) * 0.2)
Y_big <- rbind(Y1,Y1)


Ws <- myMultiCCA(xlist = list(t(X_big),t(Y_big)), penalty = c(L1,L2), type = 'standard',
                 ncomponents = 3, trace = TRUE,standardize = TRUE)
abar1_nopen_unweighted <- as.matrix(abs(Ws$ws[[1]]))
abar2_nopen_unweighted <- as.matrix(abs(Ws$ws[[2]]))
abar3_nopen_unweighted <- as.matrix(abs(Ws$ws[[3]]))
agg_smcca_nopen_unweighted <- apply(rbind(abar1_nopen_unweighted, abar2_nopen_unweighted, abar3_nopen_unweighted),1, function(x) sum(x))
plot(agg_smcca_nopen_unweighted, type = 'h')
membership_x1 <- rep(0,150)
membership_x2 <- rep(0,150)
membership_x3 <- rep(0,150)
for (i in 1:150)
{
  membership_x1[i] <- which.max(abar1_nopen_unweighted[i,])
  membership_x2[i] <- which.max(abar2_nopen_unweighted[i,])
  membership_x3[i] <- which.max(abar3_nopen_unweighted[i,1:2])
}
membership_x3 <- membership_x3[c(1:50,101:150)]
true_x1 <- c(rep(1,50), rep(3, 50), rep(2,50))
true_x3 <- c(rep(1,50), rep(2,50))
ml_test(factor(membership_x1,levels = c(1,2,3)), true_x1)
ml_test(factor(membership_x2,levels = c(1,2,3)), true_x1)
ml_test(membership_x3, true_x3)

(0.88 * 3 + 0.67 * 3 + 0.99 * 2)/8






Ws <- myMultiCCA(xlist = list(t(X1),t(X2),t(X3),t(Y1)), penalty = c(L1,L2,L3,L4), type = 'standard',
                 ncomponents = 3, trace = TRUE,standardize = TRUE, CCcoef = c(1,1,1,1,1,1))
abar1_nopen_unweighted <- as.matrix(abs(Ws$ws[[1]]))
abar2_nopen_unweighted <- as.matrix(abs(Ws$ws[[2]]))
abar3_nopen_unweighted <- as.matrix(abs(Ws$ws[[3]]))
agg_smcca_nopen_unweighted <- apply(rbind(abar1_nopen_unweighted, abar2_nopen_unweighted, abar3_nopen_unweighted),1, function(x) sum(x))
plot(agg_smcca_nopen_unweighted, type = 'h')





####################### unweighted smcca no penalty
L1 <- max(1, sqrt(nrow(X1))*0.2)
L2 <- max(1, sqrt(nrow(X2))*0.2)
L3 <- max(1, sqrt(nrow(X3))*0.2)
L4 <- max(1, sqrt(nrow(Y1))*0.2)

Ws <- myMultiCCA(xlist = list(t(X1),t(X2),t(X3),t(Y1)), penalty = c(L1,L2,L3,L4), type = 'standard',
                 ncomponents = 3, trace = TRUE,standardize = TRUE, CCcoef = c(1,1,10,1,10,10))
abar1_nopen <- as.matrix(abs(Ws$ws[[1]]))
abar2_nopen <- as.matrix(abs(Ws$ws[[2]]))
abar3_nopen <- as.matrix(abs(Ws$ws[[3]]))
agg_smcca_nopen <- apply(rbind(abar1_nopen_unweighted, abar2_nopen_unweighted, abar3_nopen_unweighted),1, function(x) sum(x))


membership_x1 <- rep(0,150)
membership_x2 <- rep(0,150)
membership_x3 <- rep(0,150)
for (i in 1:150)
{
  membership_x1[i] <- which.max(abar1_nopen[i,])
  membership_x2[i] <- which.max(abar2_nopen[i,])
  membership_x3[i] <- which.max(abar3_nopen[i,1:2])
}
membership_x3 <- membership_x3[c(1:50,101:150)]
true_x1 <- c(rep(1,50), rep(3, 50), rep(2,50))
true_x3 <- c(rep(1,50), rep(2,50))
ml_test(factor(membership_x1,levels = c(1,2,3)), true_x1)
ml_test(factor(membership_x2,levels = c(1,2,3)), true_x1)
ml_test(membership_x3, true_x3)

(0.47 * 3 + 0.45 * 3 + 0.74 * 2)/8

agg_smcca_nopen <- apply(abar1_nopen,1, function(x) sum(x))
plot(abar1_nopen[,3]~a,type = 'h')
plot(agg_smcca_nopen~a,type = 'h')
a <- list(X1,X2,X3,Y1)
save(a, file = 'experiment_data.Rdata')
