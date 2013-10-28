########################
#Q2 part d MCMC routine#
########################

library(MASS)
library(coda)

########################################################################################
########################################################################################
## Handle batch job arguments:

# 1-indexed version is used now.
args <- commandArgs(TRUE)

cat(paste0("Command-line arguments:\n"))
print(args)

####
# sim_start ==> Lowest simulation number to be analyzed by this particular batch job
###

#######################
sim_start <- 1000
length.datasets <- 200
#######################

if (length(args)==0){
  sinkit <- FALSE
  sim_num <- sim_start + 1
  set.seed(1330931)
} else {
  # Sink output to file?
  sinkit <- TRUE
  # Decide on the job number, usually start at 1000:
  sim_num <- sim_start + as.numeric(args[1])
  # Set a different random seed for every job number!!!
  set.seed(762*sim_num + 1330931)
}

# Simulation datasets numbered 1001-1200

########################################################################################
########################################################################################

data<-read.table(file=paste("data/blr_data_",as.character(sim_num),".csv", sep=""), sep=",", header=T)
data<-read.table("blr_data_1026.csv", sep=",", header=T)
m<-data$n
y<-data$y
x<-cbind(data$X1, data$X2)
beta.0<-c(0,0)
Sigma.0.inv<-diag(2)

"log.target.density" <- function(m,y,x,beta,beta.0,Sigma.0.inv){
  log.p <- t(y)%*%x%*%beta-t(m)%*%log(1+exp(x%*%beta)) - t(beta-beta.0)%*%Sigma.0.inv%*%(beta-beta.0)/2
  return(log.p)
}

"bayres.logreg"<-function(m,y,x,beta.0,Sigma.0.inv,
                          niter=10000,burnin=1000,
                          retune=100,verbose=TRUE) {

  ## Tune the covariance matrix of the proposal distribution 
  
  # Proposal distribution: beta.prop ~ N(beta.curr, V)
  # Set the initial value of V as the identity matrix
  V <- diag(2)
  beta.curr <- c(0,0)
  accept.rate=0

  while (accept.rate <= 0.3 | accept.rate >= 0.6) {
    n.accept <- 0
    for (j in 1:retune) {
      beta.prop <- mvrnorm(n=1,beta.curr,V)
      log.alpha <- log.target.density(m,y,x,beta=beta.prop,beta.0,Sigma.0.inv) - log.target.density(m,y,x,beta=beta.curr,beta.0,Sigma.0.inv)
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

  beta.curr <- c(0,0)
  # Store the samples
  beta.draws<-matrix(rep(NA, 2*(niter+burnin)), ncol=2)
  # Track tha acceptance rate:
  n.accept <- rep(0,(niter+burnin))
  
  for (i in 1:(niter+burnin)){
    beta.prop <- mvrnorm(n=1,beta.curr,V)
    log.alpha <- log.target.density(m,y,x,beta=beta.prop,beta.0,Sigma.0.inv) - log.target.density(m,y,x,beta=beta.curr,beta.0,Sigma.0.inv)
    log.u <- log(runif(1))
    if (log.u < log.alpha){
      beta.curr <- beta.prop
      n.accept[i] <- 1
    } else {
    }
    # Store the current state:
    beta.draws[i,] <- beta.curr
  }
  # Throw away burnin period
  beta.samples<-beta.draws[(burnin+1):(niter+burnin),]
  # Report the acceptance rate:
  cat(paste("Acceptance rate was ",100*round(sum(n.accept[(burnin+1):(niter+burnin)])/niter,2),"%\n",sep=""))
  
  ## Check the samples I generated
  beta1<-mcmc(beta.samples[,1])
  plot(beta1)
  beta2<-mcmc(beta.samples[,2])
  plot(beta2)
  cat(paste("ESS for beta1 = ",round(effectiveSize(beta1),4), sep=""),"\n")
  cat(paste("ESS for beta2 = ",round(effectiveSize(beta2),4), sep=""),"\n")
  cat(paste("The 95% HPD interval for beta1 = (",round(HPDinterval(beta1,0.95)[1,1],4),",", round(HPDinterval(beta1,0.95)[1,2],4),")", sep=""),"\n") 
  cat(paste("The 95% HPD interval for beta2 = (",round(HPDinterval(beta2,0.95)[1,1],4),",", round(HPDinterval(beta2,0.95)[1,2],4),")", sep=""),"\n") 

  ## Output the percentiles to a csv file 
  percentile<-matrix(rep(0, 2*99), ncol=2)
  for (i in 1:99) {
    percentile[i,1]=quantile(beta.samples[,1],probs=i/100)
    percentile[i,2]=quantile(beta.samples[,2],probs=i/100)
  }
  write.table(percentile, file=paste("results/blr_res_",as.character(sim_num),".csv", sep=""), sep=",", row.names=F, col.names=F)
  
}

bayres.logreg(m,y,x,beta.0,Sigma.0.inv,
              niter=10000,burnin=1000,
              retune=100,verbose=TRUE)