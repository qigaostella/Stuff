library(RCUDA)

cat("Loading module...\n")
m = loadModule("truncnorm.ptx")
cat("done. Extracting kernels...\n")
k_rtruncnorm = m$rtruncnorm_kernel

verbose <- TRUE
#gpu_time<-rep(0,7)
#for (i in 1:7) {
cat("done. Setting up miscellaneous stuff...\n")
N <- 10000L
rng_a <- 1L
rng_b <- 2L
rng_c <- 3L

# Truncated normal parameters:
mu <-  rep(2, N)
sigma <-  rep(1, N)
a <- rep(0, N)
b <- rep(1.5, N)

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

cat("Grid size:\n")
print(grid_dims)
cat("Block size:\n")
print(block_dims)

nthreads <- prod(grid_dims)*prod(block_dims) 
cat("Total number of threads to launch = ",nthreads,"\n")
if (nthreads < N){
  stop("Grid is not large enough...!")
}

cat("Allocating memory on device for curandStates...\n")
cu_rng_alloc_time <- system.time({
  rng_states <- cudaMalloc(elType="curandState", numEls=N, sizeof=48L) 
})

# Automated copying
cat("Launching truncated normal CUDA kernel...\n")
cu_automated_time <- system.time({
  x_temp <- rep(0,N)
  cu_rtruncnorm <- .cuda(k_rtruncnorm, "x"=x_temp, N, mu, sigma, 
                         a, b, rng_a, rng_b, rng_c, maxtries=1000L,
                         gridDim=grid_dims, blockDim=block_dims, outputs="x")
})
cat("First few values...\n")
print(head(cu_rtruncnorm))
cat(sprintf("mu = %g, sigma=%g, n=%d\n",mu[1],sigma[1],N))
cat(sprintf("Mean = %g\n",mean(cu_rtruncnorm)))
print(cu_automated_time)
#if (verbose){
#  cat(paste("RCUDA time for n = ", 10^i, ":\n", sep=""))
#  print(cu_automated_time)
#}
#gpu_time[i]<-cu_automated_time[3]
#}
#plot(gpu_time, xlab="k", main="Runtime for RCUDA code")
