rm(list=ls())
require(grid)
library(ggplot2)
library(reshape2)
library(imputation)
opt <- list(
  rate=c(1, 5, 10, 15, 20), # missing rate
  method=c("zero_impute", "rowaverage_impute", "simplyfied_knn_impute", "kNNImpute"), # impute method
  outfmt="hw3.png" #output file format
)

guess <- list(
  zero_impute=function(matrix){
    matrix[is.na(matrix)] <- 0
    matrix
  },
  
  rowaverage_impute=function(matrix){
    for(i in which(rowSums(is.na(matrix))!=0)){ 
      matrix[i, is.na(matrix[i, ])] <- mean(matrix[i, ], na.rm=T) 
    }
    matrix
  },
  
  simplyfied_knn_impute=function(matrix, k=17) {
    
    #classfied
    nomissing <- matrix[which(rowSums(is.na(matrix)) == 0), ]
    missing <- matrix[which(rowSums(is.na(matrix)) != 0), ]
    
    #euclid distance
    distance <- function(x, y) {
      common <- !is.na(x)&!is.na(y)
      sqrt(sum((x[common]-y[common])^2))
    }
    
    #impute
    for(i in 1:nrow(missing)){
      dist <- NULL
      for(j in 1:nrow(nomissing)){
        dist <- c(dist, distance(missing[i,], nomissing[j,]))
      }
      
      ksd <- order(dist)[1:k] #k_smallest distance
      
      for(id in which(is.na(missing[i,]))) {
        missing[i, id] <- sum(1/dist[ksd]*nomissing[ksd, id])/sum(1/dist[ksd])
      }
    }
    
    matrix[which(rowSums(is.na(matrix)) != 0), ] <- missing
    
    matrix
  },
  
  kNNImpute=function(matrix, k=17){
    matrix <- kNNImpute(matrix, k)$x
    
    matrix
  }
)

nrmse <- function(missing, guess, answer) { 
  i = is.na(missing)
  sqrt(mean((guess[i]-answer[i])^2)/var(answer[i])) 
}

#---

result <- matrix(0, nrow=length(opt$rate), ncol=length(opt$method), dimname=list(opt$rate, opt$method))
answer = as.matrix(read.table("ans.txt", sep="\t"))

for(i in opt$rate) {
  missing <- as.matrix(read.table(paste("m_", i, "_1.txt", sep="")), sep="\t")
  for(j in opt$method) {
    result[toString(i),j] <- nrmse(missing, guess[[j]](missing), answer)
  }
}
rm(missing, i, j)

result

#---

pic <- ggplot(melt(result, varnames=c('RATE', 'METHOD'), value.name='NRMSE'), aes(x=RATE, y=NRMSE, colour=METHOD)) + 
  geom_line(aes(group=METHOD),size=1.5) + 
  geom_point(aes(shape=METHOD),size=7) +
  xlab("Missing Rate (%)") +
  theme(axis.text=element_text(size=26),
        axis.title=element_text(size=30),
        legend.title=element_text(size=20),
        legend.text=element_text(size=16),
        legend.key.size=unit(1,"cm"),
        legend.background=element_rect(colour="gray"),
        legend.margin=unit(2.5,"cm"),
        legend.position=c(.8,.1))
pic
png(filename=opt$outfmt, width=700, height=700, bg='white')
print(pic)
dev.off()
