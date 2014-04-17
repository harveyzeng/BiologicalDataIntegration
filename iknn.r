iknn <- function(xmiss, Niter, K) {
  
  rowMeanSubstitution <- function(xmiss) {
    rowM <- rowMeans(xmiss, na.rm = T)
    rowM <- matrix(rowM, nc = ncol(xmiss), nr = nrow(xmiss), byrow=F)
    
    E <- xmiss
    E[is.na(xmiss)] <- rowM[is.na(xmiss)]
    
    return(E)
  }
  
  Eold <- rowMeanSubstitution(xmiss);
  
  miss.gene <- is.na(xmiss)
  miss.row <- which(rowSums(miss.gene)!=0)
  miss.exp <- lapply(miss.row, function(i) which(is.na(xmiss[i, ])))
 
  for(h in 1:Niter) {
    xcomplete <- Eold
    #for each target gene
    for(j in 1:length(miss.row)){
      exp <- miss.exp[[j]]
      gene <- xmiss[miss.row[j], -exp] #target gene
      cand_x <- xcomplete[-miss.row[j], -exp] #candidate matrix
      
      d <- sqrt(rowSums((cand_x - matrix(gene, nc=length(gene), nr=nrow(cand_x), byrow=T))^2))
      d.idx <- order(d)[1:K]
      d.sel <- d[d.idx]
      
      d.originId <- c(1:nrow(xcomplete))[-miss.row[j]] 
      d.originId <- d.originId[d.idx]
      y <- xcomplete[d.originId, , drop=F]
      
      const <- sum(1/d.sel)
      w <- 1/(const*d.sel)
      
      
      Eold[miss.row[j], exp] <- w%*%y[, exp, drop=F] 
    }
  }
  
  
  
  
 return(Eold) 
}

nrmse <- function(missing, guess, answer) { 
  x = is.na(missing)
  sqrt(mean((guess[x]-answer[x])^2)/var(answer[x])) 
}
answer = as.matrix(read.table("/Users/andyisman/Documents/BioInfo/lymphoma/ans.txt", sep="\t"))
missing <- as.matrix(read.table("/Users/andyisman/Documents/BioInfo/lymphoma/m_5_1.txt", sep=""), sep="\t")

cat(nrmse(missing, iknn(missing, 5, 15), answer),"\n")