##################################
#Q3 breast cancer dataset fitting#
##################################

library(MASS)
library(coda)

"log.target.density" <- function(y,x,beta,beta.0,Sigma.0.inv){
  log.p <- t(y)%*%x%*%beta-sum(log(1+exp(x%*%beta))) - t(beta-beta.0)%*%Sigma.0.inv%*%(beta-beta.0)/2
  return(log.p)
}

"logit.inv" <- function(x){
  v<-exp(x)/(1+exp(x))
  return(v)
}

data<-read.table("breast_cancer.txt", header=T)
n<-length(data$diagnosis)
p<-11
intercept<-rep(1,n)
response<-rep(0,n)
data2<-cbind(intercept, data, response)
data2[which(data2[,12]=="M"),13]<-1
cancer<-data2[,-12]

# Set the initial value of V as the covariance estimates from glm fit
fit<-glm(response~area+compactness+concavepts+concavity+fracdim+
           perimeter+radius+smoothness+symmetry+texture, family=binomial(),data=cancer)
V<-as.matrix(summary(fit)$cov.unscaled)

y<-as.matrix(cancer$response)
x<-as.matrix(cancer[1:p])

beta.0<-rep(0,p)
Sigma.0.inv<-diag(p)*0.001

"bayres.logreg"<-function(y,x,beta.0,Sigma.0.inv,
                          niter=10000,burnin=5000,
                          retune=200,verbose=TRUE) {
  
  ## Tune the covariance matrix of the proposal distribution 
  
  # Proposal distribution: beta.prop ~ N(beta.curr, V)
  beta.curr <- rep(0,p)
  accept.rate=0
  
  while (accept.rate <= 0.3 | accept.rate >= 0.6) {
    n.accept <- 0
    for (j in 1:retune) {
      beta.prop <- mvrnorm(n=1,beta.curr,V)
      log.alpha <- log.target.density(y,x,beta=beta.prop,beta.0,Sigma.0.inv) - log.target.density(y,x,beta=beta.curr,beta.0,Sigma.0.inv)
      log.u <- log(runif(1))
      if (log.u < log.alpha){
        beta.curr <- beta.prop
        n.accept <- n.accept+1
      } else {
      }
    }
    accept.rate<-round(n.accept/retune,2)
    if (accept.rate<=0.3) {
      V<-V/exp(1)
    } else 
      if (accept.rate>=0.6) {
        V<-V*exp(1)
      } else {
      }
  }
  
  cat(paste("Acceptance rate after tuning V was ",100*accept.rate,"%\n",sep=""))
  cat(paste("The covariance matrix after tuning was \n", sep=""))
  print(V)
  
  ## Begin the sampling procedure
  
  # Store the samples
  beta.draws<-matrix(rep(NA, p*(niter+burnin)), ncol=p)
  # Track tha acceptance rate:
  n.accept <- rep(0,(niter+burnin))
  
  for (i in 1:(niter+burnin)){
    beta.prop <- mvrnorm(n=1,beta.curr,V)
    log.alpha <- log.target.density(y,x,beta=beta.prop,beta.0,Sigma.0.inv) - log.target.density(y,x,beta=beta.curr,beta.0,Sigma.0.inv)
    log.u <- log(runif(1))
    if (log.u < log.alpha){
      beta.curr <- beta.prop
      n.accept[i] <- 1
    } else {
    }
    beta.draws[i,] <- beta.curr
  }
  # Throw away burnin period
  beta.samples<-beta.draws[(burnin+1):(niter+burnin),]
  # Report the acceptance rate:
  cat(paste("Acceptance rate was ",100*round(sum(n.accept[(burnin+1):(niter+burnin)])/niter,2),"%\n",sep=""))  
 
  return(beta.samples)
}  

beta.samples<-bayres.logreg(y,x,beta.0,Sigma.0.inv,
                            niter=10000,burnin=5000,
                            retune=200,verbose=TRUE)

## Check if the markov chain converges

for (j in 1:p) {
  plot(mcmc(beta.samples[,j]))
  cat(paste("ESS for ",j,"-th beta = ",round(effectiveSize(mcmc(beta.samples[,j])),4), sep=""),"\n")
}

## Compute the lag-1 autocorrelation for each component of beta
for (j in 1:p) {
  cat(paste("lag-1 autocorrelation for ",j,"-th beta = ",round(acf(beta.samples[,j],type="correlation",plot=F)$acf[2] , 4), sep=""),"\n")
}

## Create posterior credible intervals
for (j in 1:p) {
  cat(paste("The 95% posterior credible interval for ",j,"-th beta = (",round(HPDinterval(mcmc(beta.samples[,j]),0.95)[1,1],4),",", round(HPDinterval(mcmc(beta.samples[,j]),0.95)[1,2],4),")", sep=""),"\n") 
}

## Posterior predictive checking
n.new<-1000
y.new<-matrix(rep(NA, n.new*n), ncol=n)
index<-sample(1:niter, n.new)
for (i in 1:n.new) {
  beta<-beta.samples[index[i], ]
  for (j in 1:n) {
    y.new[i,j]<-rbinom(1, 1, logit.inv(x[j,]%*%beta))
  }
}
y.mean<-apply(y.new, 1, mean)
hist(y.mean, xlab="mean of the simulated datasets")
abline(v=mean(y),col="red")
legend("topright","true mean",lty=1,col="red")