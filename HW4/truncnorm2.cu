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
	 // implement the rejection sampling algorithm presented in Robert's paper
        
        if (idx < n) {
          int accept = 0;
          int numtries = 0;
          if (!isfinite(b[idx])) {
            float c = (a[idx]-mu[idx])/sigma[idx];
            float alpha = (c+sqrt(c*c+4))/2;
            while(!accept && numtries < maxtries) 
            {
              numtries = numtries+1;
              float z = c-log(curand_uniform(&rng))/alpha;
              float phi_z = exp(-(alpha-z)*(alpha-z)/2);
              float temp = mu[idx] +sigma[idx]*z;
            
              if (curand_uniform(&rng) < phi_z) {
                accept = 1;
                x[idx] = temp;
              }
            }
          } else {
            float c = -(b[idx]-mu[idx])/sigma[idx];
            float alpha = (c+sqrt(c*c+4))/2;
            while(!accept && numtries < maxtries) 
            {
              numtries = numtries+1;
              float z = c-log(curand_uniform(&rng))/alpha;
              float phi_z = exp(-(alpha-z)*(alpha-z)/2);
              float temp = mu[idx] - sigma[idx]*z;
            
              if (curand_uniform(&rng) < phi_z) {
                accept = 1;
                x[idx] = temp;
              }
          }
  
            if (numtries >= maxtries) {
              printf("Reach Maximum Number of Tries. \n");
            }
          }
        }
        
        return;
    }

    } // END extern "C"