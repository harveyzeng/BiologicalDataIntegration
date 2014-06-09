
IKNN <- function(m, type="d", k=15, iteration=2){
  sim <- list(
    d = function(a,b){
      1/sqrt(sum((a-b)**2))
    },
    r = function(a,b){
      abs(cor(a,b,method="pearson"))		 
    },
    t = function(a,b){
      abs(vcosw(a,b))	
    }
  )
  complete <- m
  incomplete <- m[which(0!=rowSums(is.na(m))),]
  na.row <- which(rowSums(is.na(m))!=0)
  complete <- t(apply(complete,1,function(y){
    y[is.na(y)] <- mean(y,na.rm=T);
    y
  }))
  for(j in seq(iteration)){ # do KNN for "iteration" times
    filled <- apply(incomplete, 1, function(y){
      na.vector <- is.na(y)
      nna <- which(!na.vector)
      na <- which(na.vector)
      predition.feature <- y[nna]
      simularity <- apply(as.matrix(complete[,nna]),1,function(training.feature){sim[[type]](training.feature,predition.feature)})
      training.answer <- as.matrix(complete[,na])
      rank <- order(-simularity)[2:(k+1)]
      simularity.k <- simularity[rank]
      predition.answer <- colSums(as.matrix(simularity.k*training.answer[rank,]))/sum(simularity.k)
      y[na] <- predition.answer
      y
    })
    complete[na.row,] <- t(filled) # update the complete 
  }
  complete
}
IKNNimpute<-function(x, k=15, sim.method="EuDist", iter=2, e=1e-3) {
  similarityCal<-function(vec, mat, method="EuDist"){
    methods<-c("EuDist","cor","Angle")
    switch(match.arg(method,methods),
           EuDist=1/sqrt(rowSums((mat-matrix(vec,nc=length(vec),nr=nrow(mat),byrow=TRUE))^2)),
           cor=apply(mat,1,function(i) abs(cor(i,vec,use="everything",method="pearson"))),
           Angle=apply(mat,1,function(i) abs(sum(i * vec)/sqrt(sum(i^2)*sum(vec^2))))
    )
  }
  
  cat(sprintf("k = %g, iter= %g, e = %g\n",k,iter,e))
  missIdx<-is.na(x)
  rowNum<-nrow(x)
  miss.RowIdx<-which(rowSums(missIdx) != 0)
  RAVGimpute <- function(xmiss) {
    rowM <- rowMeans(xmiss, na.rm = T)
    rowM <- matrix(rowM, nc = ncol(xmiss), nr = nrow(xmiss), byrow=F)
    
    E <- xmiss
    E[is.na(xmiss)] <- rowM[is.na(xmiss)]
    
    return(E)
  }
  x.ravged<-RAVGimpute(x)
  cat("Row average imputation completed!\n")
  x.miss<-(cbind(1:rowNum, x))[miss.RowIdx,]
  x.complete<-x.ravged
  cat("The size of x.complete:",dim(x.complete),"\n")
  err<-99
  for (r in 1:iter) {
    x.old<-x.complete
    err.old<-err
    cat(sprintf("Start the %g cycle of iknn imputation\n", r))
    x.imputed<-t(apply(x.miss, 1, function(j) {
      rowIdx<-j[1]
      j.origin<-j[-1]
      neighbor.pool<-x.complete[-rowIdx,]
      target<-x.complete[rowIdx,]
      dist.list<-similarityCal(target, neighbor.pool, method=sim.method)
      neighborsIdx<-order(dist.list,decreasing=T)[1:k]
      missColIdx<-which(is.na(j.origin))
      estimation<-sapply(missColIdx, function(h){
        weight<-dist.list[neighborsIdx]
        weightedAvg<-weight %*% neighbor.pool[neighborsIdx,h]/sum(weight)
        return(weightedAvg)
      })
      j.origin[missColIdx]<-estimation
      return(j.origin)
    }))
    cat(sprintf("The %g cycle imputation has been completed\n",r))
    x.complete[miss.RowIdx,]<-x.imputed
    err<-sum((x.old[missIdx]-x.complete[missIdx])^2)
    cat("err:",err,"\n")
    if (r>1) if ((err.old/err < 4) | (err < e)) break
  }
  x<-x.complete
  return(x)
}

iknn_1 <- function(xmiss, K=15, sim.method="EuDist", Niter=3){
  similarityCal<-function(vec, mat, method="EuDist"){
    methods<-c("EuDist","cor","Angle")
    switch(match.arg(method,methods),
           EuDist=1/sqrt(rowSums((mat-matrix(vec,nc=length(vec),nr=nrow(mat),byrow=TRUE))^2)),
           cor=apply(mat,1,function(i) abs(cor(i,vec,use="everything",method="pearson"))),
           Angle=apply(mat,1,function(i) abs(sum(i * vec)/sqrt(sum(i^2)*sum(vec^2))))
    )
  }
  
  rowMeanSubstitution <- function(xmiss) {
    rowM <- rowMeans(xmiss, na.rm = T)
    rowM <- matrix(rowM, nc = ncol(xmiss), nr = nrow(xmiss), byrow=F)
    
    E <- xmiss
    E[is.na(xmiss)] <- rowM[is.na(xmiss)]
    
    return(E)
  }
  
  xcomplete <- rowMeanSubstitution(xmiss);
  
  miss.gene <- is.na(xmiss)
  miss.row <- which(rowSums(miss.gene)!=0)
  miss.exp <- lapply(miss.row, function(i) which(is.na(xmiss[i, ])))
  
  for(h in 1:Niter) {
    xcomplete[miss.row, ] <- t(apply(xmiss[miss.row, ], 1, function(row){
      row.miss <- is.na(row)
      row.exp <- which(row.miss)
      d <- similarityCal(row[-row.exp], xcomplete[, -row.exp, drop=F], sim.method)
      id.idx <- order(d, decreasing=T)[2:(K+1)]
      id.sel <- d[id.idx]
      const <- sum(1/id.sel)
      w <- 1/(const*id.sel)
      row[row.exp] <- w %*% xcomplete[id.idx, row.exp, drop=F]
      return (row)
    }))
  }
  return(xcomplete) 
}


nrmse <- function(missing, guess, answer) { 
  x = is.na(missing)
  sqrt(mean((guess[x]-answer[x])^2)/var(answer[x])) 
}
answer = as.matrix(read.table("/Users/andyisman/Documents/BioInfo/lymphoma/ans.txt", sep="\t"))
missing <- as.matrix(read.table("/Users/andyisman/Documents/BioInfo/lymphoma/m_5_1.txt", sep=""), sep="\t")

ptm <- proc.time()
cat(nrmse(missing, iknn(missing), answer),"\n", proc.time() - ptm, "\n")

