    #include <stdio.h>
    #include <stdlib.h>
    #include <cuda.h>
    #include <curand_kernel.h>
    #include <math_constants.h>

    extern "C"
    {

    __global__ void
    rtruncnorm_kernel(
      float *x,      // Vector to contain returned samples 
      int n,         // Number of samples to return
      float *mu,     // Vector of mu's
      float *sigma,  // Vector of sigma's
      float *a,      // Vector of lower-truncation values
      float *b,      // Vector of upper-truncation values
      int rng_a,     // RNG seed constant
      int rng_b,     // RNG seed constant
      int rng_c,      // RNG seed constant
      int maxtries)     
    {
        // Usual block/thread indexing...
        int myblock = blockIdx.x + blockIdx.y * gridDim.x;
        int blocksize = blockDim.x * blockDim.y * blockDim.z;
        int subthread = threadIdx.z*(blockDim.x * blockDim.y) + threadIdx.y*blockDim.x + threadIdx.x;
        int idx = myblock * blocksize + subthread;
        
        // Initialize RNG 
        curandState rng; // random number generator
        curand_init(rng_a+idx*rng_b, rng_c, 0, &rng); // rng_c sequence number (related to iteration number in MCMC)
        // initialize different seed for different index
        
        // Sample the truncated normal:
        // mu for this index is mu[idx]
        // sigma for this index is sigma[idx]
        // a for this index is a[idx]
        // b for this index is b[idx]
        
        // X_i ~ Truncated-Normal(mu_i, sigma_i; [a_i, b_i])
        
        // Sample N(mu, sigma^2):
        if (idx < n) {
          int accept = 0;
          int numtries = 0;
          while(!accept && numtries < maxtries) 
          {
            float temp = mu[idx] +sigma[idx]*curand_normal(&rng);
            numtries = numtries+1;
            if (temp>=a[idx] && temp<=b[idx]) {
              accept = 1;
              x[idx] = temp;
            }
          }
  
          if (numtries >= maxtries) {
           printf("Reach Maximum Number of Tries. \n");
          }
        }

        return;
    }

    } // END extern "C"