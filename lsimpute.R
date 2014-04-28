LsImpute <- function(xmiss, K){
  
  miss.gene <- is.na(xmiss)
  miss.row <- which(rowSums(miss.gene)!=0)
  
  xcomplete <- xmiss[-miss.row, ]
  xincomplete <- xmiss[miss.row, ]

  impute <- function(row) {
    
    row.exp <- which(is.na(row))
    gene <- row[-row.exp]
    cand_x <- xcomplete[, -row.exp]
    
    sim <- apply(cand_x, 1, function(x) {abs(cor(x, gene, method="pearson"))})
    row.idx <- order(sim, decreasing=T)[1:K]
    row.r <- sim[row.idx]
    row.cand <- cand_x[row.idx, ]
    lg <- apply(row.cand, 1, function(x){lm(x~gene)$coefficients})
    row.impcand <- xcomplete[row.idx, row.exp]
  
    y <- matrix(0, nc=ncol(row.impcand), nr=nrow(row.impcand))
    w <- (row.r**2/(1-row.r**2+0.000001))**2
    sw <- sum(w)
    w <- w/sw

    for(i in 1:nrow(row.impcand)) {
      y[i, ] <- lg[2, i] * row.impcand[i, ] + lg[1, i]
    }
    row[row.exp]<-apply(y, 2, function(x){sum(w*x)})
    return(row)
  }
  
  xmiss[miss.row, ] <- t(apply(xincomplete, 1, impute))
 
  return (xmiss)
}

nrmse <- function(missing, guess, answer) { 
  x = is.na(missing)
  sqrt(mean((guess[x]-answer[x])^2)/var(answer[x])) 
}
answer = as.matrix(read.table("/Users/andyisman/Documents/BioInfo/lymphoma/ans.txt", sep="\t"))
missing <- as.matrix(read.table("/Users/andyisman/Documents/BioInfo/lymphoma/m_10_1.txt", sep=""), sep="\t")

cat(nrmse(missing, LsImpute(missing, 15), answer))