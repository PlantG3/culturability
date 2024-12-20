size.factor.DESeq <- function(counts) {
### compute size factor using the method described in the DESeq paper
### citation: Genome Biology 2010, 11:R106
### Sanzhen Liu
### 8/14/2014

###########################################
  geometric.mean <- function (x) {
  ### function to compute geometric mean
    gm <- 0
    if (sum(x==0)<1) {
      gm <- exp(mean(log(x)))
    }
    return (gm)
  }
###########################################

  gmvalue <- apply(counts, 1, geometric.mean) 
  sf <- rep(NA, ncol(counts))  ### initiate size factors
  for (i in 1:ncol(counts)) {
    cur.counts <- counts[, i]
    sf[i] <- median(cur.counts[gmvalue>0]/gmvalue[gmvalue>0])
  }
  return(sf)
}
