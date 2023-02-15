# simulate data
set.seed(1234) 

# correlation structure
index_4_way <- matrix(0, nrow = 3, ncol = 300)
index_124 <- matrix(0, nrow = 2, ncol = 300)
index_134 <- matrix(0, nrow = 2, ncol = 300)
index_234 <- matrix(0, nrow = 2, ncol = 300)
index_14 <- matrix(0, nrow = 2, ncol = 300)
index_24 <- matrix(0, nrow = 2, ncol = 300)
index_34 <- matrix(0, nrow = 2, ncol = 300)

index_123 <- matrix(0, nrow = 3, ncol = 300)
index_12 <- matrix(0, nrow = 2, ncol = 300)
index_13 <- matrix(0, nrow = 2, ncol = 300)
index_23 <- matrix(0, nrow = 2, ncol = 300)

index_4_way[,1:20] <- runif(20*3,0.4, 0.6) * 2 * (rbinom(20*3, 1, 0.5) - 0.5)
index_124[,21:40] <- runif(20*2,0.4, 0.6) * 2 * (rbinom(20*2, 1, 0.5) - 0.5)
index_134[,41:60] <- runif(20*2,0.4, 0.6) * 2 * (rbinom(20*2, 1, 0.5) - 0.5)
index_234[,61:80] <- runif(20*2,0.4, 0.6) * 2 * (rbinom(20*2, 1, 0.5) - 0.5)
index_14[,81:100] <- runif(20*2,0.4, 0.6) * 2 * (rbinom(20*2, 1, 0.5) - 0.5)
index_24[,81:100] <- runif(20*2,0.4, 0.6) * 2 * (rbinom(20*2, 1, 0.5) - 0.5)
index_34[,81:100] <- runif(20*2,0.4, 0.6) * 2 * (rbinom(20*2, 1, 0.5) - 0.5)

# define number of observations
n_obs <- 100
# generate latent variables
latent <- MASS::mvrnorm(n = n_obs, mu = rep(0, 8), Sigma = diag(8))

# simulate data matrix
x1 <- latent[,1] %*% t(index_4_way[1,]) + 
  latent[,2] %*% t(index_124[1,]) + latent[,3] %*% t(index_134[1,]) +
  latent[,5] %*% t(index_14[1,]) + 0.2 * matrix(rnorm(300*n_obs,0,1),nrow = n_obs, ncol = 300)

x2 <- latent[,1] %*% t(index_4_way[2,]) + 
  latent[,2] %*% t(index_124[2,]) + latent[,4] %*% t(index_234[1,]) +
  latent[,6] %*% t(index_24[1,]) + 0.2 * matrix(rnorm(300*n_obs,0,1),nrow = n_obs, ncol = 300)


x3 <- latent[,1] %*% t(index_4_way[3,]) + 
  latent[,3] %*% t(index_134[2,]) + latent[,4] %*% t(index_234[2,]) +
  latent[,7] %*% t(index_34[1,]) + 0.2 * matrix(rnorm(300*n_obs,0,1),nrow = n_obs, ncol = 300)

# scale data
x4 <- scale(rowSums(latent[,1:7])/7 + 0.2 * latent[,8])
x1 <- scale(x1)
x2 <- scale(x2)
x3 <- scale(x3)
# define feature names for simulated data
colnames(x1) <- paste0('gene_',1:300)
colnames(x2) <- paste0('protein_',1:300)
colnames(x3) <- paste0('metabolite_',1:300)

# compile all multi-omics data into data list
data <- list(x1, x2, x3, x4)
# define data type for data list
dtype <- c('gene', 'protein', 'metabolite', 'pheno')

# run SGTCCA-Net
result <- SGTCCA_Net(data_list = data, correlation_list = list(c(1,2,3,4), c(1,2,4)), network_exclude_data = 4, data_type = dtype,
           pheno = x4, max_mod_size = 100, saving_dir = 'C:/Users/liux4/Documents')

# save global network result
save(result, file = 'C:/Users/liux4/Documents/gtcca_result.Rdata') 