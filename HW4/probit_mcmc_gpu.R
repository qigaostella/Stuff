#library(mvtnorm)
#library(truncnorm)
#library(coda)
library(MASS)

# Define a Function Which Calculates Inverve Using SVD
inv<-function(X) {
  u<-svd(X)$u
  d<-svd(X)$d
  v<-svd(X)$v
  inv<-v%*%diag(1/d)%*%t(u)
  return(inv)
}

library(RCUDA)
cat("Loading module...\n")
m = loadModule("truncnorm2.ptx")
cat("done. Extracting kernels...\n")
k_rtruncnorm = m$rtruncnorm_kernel
rng_a <- 8L
rng_b <- 26L
data<-read.table("mini_data.txt", sep="", header=T)
y<-as.matrix(data[,1])
X<-as.matrix(data[,2:9])
colnames(X)<-NULL
N<-as.integer(dim(X)[1])
p<-dim(X)[2]
beta_0<-rep(0.0, p)
Sigma_0_inv<-matrix(rep(0.0, p*p), ncol=p)

"compute_grid" <- function(N,sqrt_threads_per_block=16L,grid_nd=1)
{
  block_dims <- c(as.integer(sqrt_threads_per_block), as.integer(sqrt_threads_per_block), 1L)
  threads_per_block <- prod(block_dims)
  if (grid_nd==1){
    grid_d1 <- as.integer(max(1L,ceiling(N/threads_per_block)))
    grid_d2 <- 1L
  } else {
    grid_d1 <- as.integer(max(1L, floor(sqrt(N/threads_per_block))))
    grid_d2 <- as.integer(ceiling(N/(grid_d1*threads_per_block)))
  }
  grid_dims <- c(grid_d1, grid_d2, 1L)
  return(list("grid_dims"=grid_dims,"block_dims"=block_dims))
}

# Fix the grid/block dimensions
bg <- compute_grid(N)
grid_dims <- bg$grid_dims
block_dims <- bg$block_dims

nthreads <- prod(grid_dims)*prod(block_dims) 
cat("Total number of threads to launch = ",nthreads,"\n")
if (nthreads < N){
  stop("Grid is not large enough...!")
}

cat("Allocating memory on device for curandStates...\n")
cu_rng_alloc_time <- system.time({
  rng_states <- cudaMalloc(elType="curandState", numEls=N, sizeof=48L) 
})

probit_mcmc_gpu<-function(y, X, beta_0, Sigma_0_inv, 
                          niter, burnin, block_dims, grid_dims) {
  n<-dim(X)[1]
  p<-dim(X)[2]
  
  beta_draws<-matrix(rep(NA, p*(niter+burnin)), ncol=p)
  beta_draws[1,]<-beta_0
  Sigma_post<-inv(Sigma_0_inv+t(X)%*%X)
  
  sigma <-  rep(1.0, N)
  a<-rep(0.0, N)
  b<-rep(0.0, N)
  a[which(y==0)]=-Inf
  b[which(y==1)]=Inf
  
  for (i in 1:(burnin+niter-1)) {
    rng_c <- as.integer(i)
    mu <- X%*%beta_draws[i,]
    x_temp <- rep(0.0,N)
    z <- .cuda(k_rtruncnorm, "x"=x_temp, N, mu, sigma, 
               a, b, rng_a, rng_b, rng_c, maxtries=1000L,
               gridDim=grid_dims, blockDim=block_dims, outputs="x")
    beta_post<-Sigma_post%*%(Sigma_0_inv%*%beta_0+t(X)%*%z)
    beta_draws[i+1,]<-mvrnorm(1, beta_post, Sigma_post)
    if(i<=burnin && i%%100==0) cat(i," iterations within burnin completed!\n", sep="")
    if(i>burnin && i%%100==0) cat((i-burnin)," iterations after burnin completed!\n", sep="")
  }
  
  cat("MCMC sampling completed!\n", sep="")
  
  return(beta_draws[(burnin+1):(niter+burnin),])
  
}

gpu_time <- system.time({
  beta_samples<-probit_mcmc_gpu(y, X, beta_0, Sigma_0_inv, niter=5000, burnin=2000, block_dims, grid_dims)
})
cat("The runtime of gpu code for mini data is:\n")
print(gpu_time)
true_pars<-as.matrix(read.table("mini_pars.txt", sep="", header=T))
colnames(true_pars)<-NULL
cat("The true beta is:\n")
print(t(true_pars))
post_mean<-colMeans(beta_samples)
cat("The posterior mean of beta is:\n")
print(post_mean)

#Traceplots
beta_samples<-read.table("beta_mini.txt", header=F)
par(mfrow=c(4,2))
for (i in 1:p) {
  plot(beta_samples[,i], type="l", ylab="beta_samples", main=paste("Traceplot of Var",i,sep=""))
  abline(h=true_pars[i], col="red")
}

#dev.off()