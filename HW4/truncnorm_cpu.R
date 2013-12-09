library(tmvtnorm)
library(truncnorm)

n=10000
TN_samples<-rtruncnorm(n, 0, 1.5, 2, 1)
hist(TN_samples, xlab="Truncated Normal Samples", main="Histogram of CPU-generated TN(2,1;(0,1.5)) Samples")
sim_mean<-mean(TN_samples) 
theo_mean<-mtmvnorm(2, 1, 0, 1.5)$tmean # 0.9570067

TN_samples<-rtruncnorm(n, -Inf, 1.5, 2, 1)
hist(TN_samples, xlab="Truncated Normal Samples", main="Histogram of CPU-generated TN(2,1;(-Inf,1.5)) Samples")
sim_mean<-mean(TN_samples)
theo_mean<-mtmvnorm(2, 1, -Inf, 1.5)$tmean # 0.8589222

TN_samples<-rtruncnorm(n, 0, Inf, 2, 1)
hist(TN_samples, xlab="Truncated Normal Samples", main="Histogram of CPU-generated TN(2,1;(0,Inf)) Samples")
sim_mean<-mean(TN_samples) 
theo_mean<-mtmvnorm(2, 1, 0, Inf)$tmean # 2.055248

TN_samples<-rtruncnorm(n, -Inf, -10, 0, 1)
hist(TN_samples, xlab="Truncated Normal Samples", main="Histogram of CPU-generated TN(0,1;(-Inf,-10)) Samples")
sim_mean<-mean(TN_samples) 
theo_mean<-mtmvnorm(0, 1, -Inf, -10)$tmean # -10.09809

verbose <- TRUE
cpu_time<-rep(0,8)
for (i in 1:8) {
  r_time <- system.time({
    r_samples <- rtruncnorm(n=10^i, 0, 1.5, 2, 1)
  })
  if (verbose){
    cat(paste("R time for n = ", 10^i, ":\n", sep=""))
    print(r_time)
  }
  cpu_time[i]<-r_time[3]
}
plot(cpu_time, xlab="k", main="Runtime for Pure R code")