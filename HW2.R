library(reshape2)
library(missForest)
library(ggplot2)
require(grid)
#read data
ans <- read.table("/Users//andyisman//Documents/BioInfo/lymphoma/ans.txt", sep='\t')

data = list()
ndata = list()

for(i in c(1,5,10,15,20)){
  data[[length(data) + 1]] <- read.table(paste("/Users/andyisman/Documents//BioInfo//lymphoma//","m_",i,"_1.txt", sep=""), sep='\t')
  ndata[[length(ndata) + 1]] <- is.na(data[[length(data)]])
}


zeroImpute <- function( rawdata ) {
  rawdata[ is.na( rawdata ) ] <- 0;
  return ( rawdata )
}

averageImpute <- function( rawdata ) {
  impute <- function( x ){
    x[ is.na( x ) ] <- mean( x, na.rm = T )
    return ( x )
  }
  rawdata <- t( apply( rawdata, 1, impute ) )
}

#impute 
idata = list()
for(i in 1:length(data)) {
  idata[[ i ]] <- list(zeroImpute(data[[i]]), averageImpute(data[[i]]))
}


a1 <- matrix(c(1, 'zeroimpute', nrmse(idata[[i]][[1]], data[[i]],ans)), nc=3)
colnames(a1) <- c('M','Method','NRMSE')
a1 <- rbind(a1, (c(1, 'aveimpute', nrmse(idata[[i]][[2]], data[[i]],ans))))

for(i in 2:length(data)) {
  a1<- rbind(a1, c((i-1)*5, 'zeroimpute', nrmse(idata[[i]][[1]], data[[i]], ans)))
  a1<- rbind(a1, c((i-1)*5, 'aveimpute', nrmse(idata[[i]][[2]], data[[i]],ans)))
}

a2<-data.frame(M=as.integer(a1[,1]),Method=a1[,2],NRMSE=as.numeric(a1[,3]))

p1<-ggplot(a2,aes(x=M,y=NRMSE,colour=Method)) + 
  geom_line(aes(group=Method),size=1.5) + 
  geom_point(aes(shape=Method),size=7) +
  xlab("Missing Rate (%)") +
  theme(axis.text=element_text(size=26),
        axis.title=element_text(size=30),
        legend.title=element_text(size=28),
        legend.text=element_text(size=26),
        legend.key.size=unit(1,"cm"),
        legend.background=element_rect(colour="gray"),
        legend.margin=unit(2.5,"cm"),
        legend.position=c(0.8,0.2))





  