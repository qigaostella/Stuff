mini <- FALSE

#============================== Setup for running on Gauss... ==============================#

args <- commandArgs(TRUE)

cat("Command-line arguments:\n")
print(args)

####
# sim_start ==> Lowest possible dataset number
####

###################
sim_start <- 1000
###################

if (length(args)==0){
  sim_num <- sim_start + 1
  set.seed(121231)
} else {
  # SLURM can use either 0- or 1-indexing...
  # Lets use 1-indexing here...
  sim_num <- sim_start + as.numeric(args[1])
  sim_seed <- (762*(sim_num-1) + 121231)
}

cat(paste("\nAnalyzing dataset number ",sim_num,"...\n\n",sep=""))

# Find r and s indices:

s_index<-ceiling((sim_num-1000)/50)
r_index<-(sim_num-1000)%%50
if (r_index==0) r_index<-50

#============================== Run the simulation study ==============================#

# Load packages:
library(BH)
library(bigmemory.sri)
library(bigmemory)
library(biganalytics)

# I/O specifications:
datapath <- "/home/pdbaines/data"
outpath <- "output/"

# mini or full?
if (mini){
  rootfilename <- "blb_lin_reg_mini"
} else {
  rootfilename <- "blb_lin_reg_data"
}


##########################
#Bag of Little Bootstraps#
##########################

# Read in the data
descriptorfilename <- paste0(rootfilename,".desc")
descriptorfile <- paste(datapath,descriptorfilename,sep="/")
dat<-attach.big.matrix(dget(descriptorfile),backingpath=datapath)

# Find sub-sample size b
n<-nrow(dat)
gamma<-0.7
b<-floor(n^gamma)

# Sample b rows from the full dataset
index<-{set.seed(s_index); sample(1:n, b, replace=F)}
subsample<-as.data.frame(dat[index,])
# Draw a bootstrap sample of size n from the subsample
rep.times<-{set.seed(sim_num);rmultinom(1,n, prob=rep(1, b)/b)}

# Fit a weighted least squares regression
if (mini){
  fit<-lm(subsample$V41~.-1, data=subsample, weights=rep.times)
} else {
  fit<-lm(subsample$V1001~.-1, data=subsample, weights=rep.times)
}
coef<-fit$coef

# Save estimates to file
outfile=paste0("output/", "coef_", sprintf("%02d", s_index), "_", 
              sprintf("%02d", r_index), ".txt")
write.table(coef, file=outfile, sep=",", row.names=F, col.names=F)

# Produce an index plot
se<-read.table("blb_lin_reg_SE.txt", header=TRUE)
plot(1:1000, se[,1], xlab="index", ylab="SE", main="Index Plot")
abline(h=0.01, col="red")