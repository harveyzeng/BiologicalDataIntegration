
iknn <- function(xmiss, Niter, K) {
  
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
  xincomplete <- xmiss[miss.row, ]
  
  for(h in 1:Niter) {
    #for each target gene
    impute <- function(row) {
      row.miss <- is.na(row)
      row.exp <- which(row.miss)
      gene <- row[-row.exp]
      cand_x <- xcomplete[, -row.exp]
      
      d <- sqrt(rowSums((cand_x - matrix(gene, nc=length(gene), nr=nrow(cand_x), byrow=T))^2))
      id.idx <- order(d)[2:(K+1)]
      id.sel <- d[id.idx]
      
      const <- sum(1/id.sel)
      w <- 1/(const*id.sel)
      
      row[row.exp] <- w %*% xcomplete[id.idx, row.exp, drop=F]
      return (row)
    }
    xcomplete[miss.row, ] <- t(apply(xincomplete, 1, impute))
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
cat(nrmse(missing, iknn(missing, 2, 15), answer),"\n", proc.time() - ptm, "\n")

