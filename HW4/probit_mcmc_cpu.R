#library(mvtnorm)
library(truncnorm)
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

probit_mcmc_cpu<-function(y, X, beta_0,
                          Sigma_0_inv, niter, burnin) {
  n<-dim(X)[1]
  p<-dim(X)[2]
  
  z<-rep(0, n)
  l<-length(y[which(y==0)])
  a<-rep(0.0, n)
  b<-rep(0.0, n)
  a[which(y==0)]=-Inf
  b[which(y==1)]=Inf
  
  beta_draws<-matrix(rep(NA, p*(niter+burnin)), ncol=p)
  beta_draws[1,]<-beta_0
  Sigma_post<-inv(Sigma_0_inv+t(X)%*%X)
  
  for (i in 1:(burnin+niter-1)) {
    #z<-sample_z(beta_draws[i,])
    z<-rtruncnorm(n, a, b, X%*%beta_draws[i,], rep(1,n))
    beta_post<-Sigma_post%*%(Sigma_0_inv%*%beta_0+t(X)%*%z)
    beta_draws[i+1,]<-mvrnorm(1, beta_post, Sigma_post)
    if(i<=burnin && i%%100==0) cat(i," iterations within burnin completed!\n", sep="")
    if(i>burnin && i%%100==0) cat((i-burnin)," iterations after burnin completed!\n", sep="")
  }
  
  cat("MCMC sampling completed!\n", sep="")
  
  return(beta_draws[(burnin+1):(niter+burnin),])
  
}

data<-read.table("mini_data.txt", sep="", header=T)
y<-as.matrix(data[,1])
X<-as.matrix(data[,2:9])
glm.coef<-glm(y~.-1, data=data, family=binomial(link="probit"))$coefficients
print(glm.coef)
  
colnames(X)<-NULL
n<-dim(X)[1]
p<-dim(X)[2]
beta_0<-rep(0, p)
Sigma_0_inv<-matrix(rep(0, p*p), ncol=p)
cpu_time <- system.time({
beta_samples<-probit_mcmc_cpu(y, X, beta_0, Sigma_0_inv, niter=5000, burnin=2000)
})
cat("The runtime of cpu code for mini data is:\n")
print(cpu_time)
true_pars<-as.matrix(read.table("mini_pars.txt", sep="", header=T))
colnames(true_pars)<-NULL
cat("The true beta is:\n")
print(t(true_pars))
post_mean<-colMeans(beta_samples)
cat("The posterior mean of beta is:\n")
print(post_mean)

#Traceplots
par(mfrow=c(4,2))
for (i in 1:p) {
  plot(beta_samples[,i], type="l", ylab="beta_samples", main=paste("Traceplot of Var",i,sep=""))
  abline(h=true_pars[i], col="red")
}
#dev.off()