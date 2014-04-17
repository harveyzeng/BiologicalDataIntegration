sknn <- function(xmiss, K) {
  
  miss.gene <- is.na(xmiss)
  miss.row <- which(rowSums(miss.gene) != 0)
  miss.exp <- lapply(miss.row, function(i) which(is.na(xmiss[i, ])))
 
  xincomplete <- xmiss[miss.row, ]
  xcomplete <- xmiss[-miss.row, ]

  xincomplete <- xincomplete[order(rowSums(xincomplete)), ]
  
  for(i in 1:nrow(xincomplete)){
    exp <- miss.exp[[i]]
    gene <- xincomplete[i, -exp]
    cand_x <- xcomplete[, -exp]
    
    d <- sqrt(rowSums((cand_x - matrix(gene, nc = length(gene), nr = nrow(cand_x), byrow=T))^2))
    
    id.idx <- order(d)[1:K]
    id.sel <- d[id.idx]
    
    const <- sum(1/id.sel)
    w <- 1/(const*id.sel)
    
    xincomplete[i, exp] <- w %*% xcomplete[id.idx, exp, drop=F]
    xcomplete <- rbind(xcomplete, xincomplete[i, ])
    xmiss[miss.row[i], ] <- xincomplete[i, ]
  }
  
  return (xmiss)
}

nrmse <- function(missing, guess, answer) { 
  x = is.na(missing)
  sqrt(mean((guess[x]-answer[x])^2)/var(answer[x])) 
}
answer = as.matrix(read.table("/Users/andyisman/Documents/BioInfo/lymphoma/ans.txt", sep="\t"))
missing <- as.matrix(read.table("/Users/andyisman/Documents/BioInfo/lymphoma/m_5_1.txt", sep=""), sep="\t")

cat(nrmse(missing, sknn(missing, 15), answer),"\n")