



sknn <- function(xmiss, K) {
  
  miss.gene <- is.na(xmiss)
  miss.row <- which(rowSums(miss.gene) != 0)
  xincomplete <- xmiss[miss.row, ]
  miss.inc <- is.na(xincomplete)
  miss.origin <- order(rowSums(miss.inc))
  xincomplete <- xincomplete[miss.origin, ]
  xcomplete <<- xmiss[-miss.row, ]
  xtmp <- matrix(nc=ncol(xincomplete))
  impute <- function(row) {  
    row.miss <- is.na(row)
    row.exp <- which(row.miss)
    gene <- row[-row.exp]
    cand_x <- xcomplete[, -row.exp]
    d <- sqrt(rowSums((cand_x - matrix(gene, nc = length(gene), nr = nrow(cand_x), byrow=T))^2))
    
    id.idx <- order(d)[1:K]
    id.sel <- d[id.idx]
    
    const <- sum(1/id.sel)
    w <- 1/(const*id.sel)
    
    row[row.exp] <- w%*%xcomplete[id.idx, row.exp, drop=F]
    xcomplete <<- rbind(xcomplete, row)
    
    return (row)
  } 
  
  xmiss[miss.row[miss.origin], ] <- t(apply(xincomplete, 1, impute))
    
  return (xmiss)
}

nrmse <- function(missing, guess, answer) { 
  x = is.na(missing)
  sqrt(mean((guess[x]-answer[x])^2)/var(answer[x])) 
}
answer = as.matrix(read.table("/Users/andyisman/Documents/BioInfo/lymphoma/ans.txt", sep="\t"))
missing <- as.matrix(read.table("/Users/andyisman/Documents/BioInfo/lymphoma/m_5_1.txt", sep=""), sep="\t")
ptc <- proc.time()
cat(nrmse(missing, sknn(missing, 15), answer),"\n")
