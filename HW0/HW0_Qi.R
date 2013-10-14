library(utils)

# Q1
####
x<-c(1:100)

for (i in 1:floor(100/3)) {
  x[3*i]<-"Fizz"
}

for (i in 1:floor(100/5)) {
  x[5*i]<-"Buzz"
}

for (i in 1:floor(100/15)) {
  x[15*i]<-"FizzBuzz"
}

x

# Q2
####
x<-runif(10000, min = 0, max = 2*pi)
y<-runif(10000, min = 0, max = 1)
u<-y*cos(x)
v<-y*sin(x)
plot(u, v, main="Scatterplot of the (u,v) Pairs")
r<-sqrt(u^2+v^2)
plot(r, main="Plot of r")
# r is a uniform distribution between 0 and 1.

# Q3
####
s<-"Hello, my name is Bob. I am a statistician. I like statistics very much."
s.split<-unlist(strsplit(s, split=NULL))
l=length(s.split)
for (i in 1:l) {
  filename = paste("out_", sprintf("%02d", i), ".txt", sep="")
  write.table(s.split[i], filename, quote=T, row.names = F, col.names = F)
}

file_list<-list.files(pattern = "^out")

dataset<-read.table(file_list[1], header=F, sep="")

for (j in 2:l){
  #temp_dataset <-read.table(file_list[j], header=F, sep="", blank.lines.skip = F)
  temp_dataset<-scan(file_list[j],what="character", quiet=TRUE)
  ##Scan is a better function than read.table in this case!
  dataset<-cbind(dataset, temp_dataset)
  #write.table(temp_dataset, "output.txt", append = T, quote=F, sep="", 
  #            na = " ", row.names = F, col.names = F)
  rm(temp_dataset)
}

write.table(dataset, "sentence.txt",quote=F, row.names = F, col.names = F, sep="")

# Q6
####
