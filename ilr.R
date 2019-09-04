ilr <- function(ID, M){
  library(MASS)
  n <- dim(M)[1]
  p <- dim(M)[2]
  M1 <- t(t(M[,1]))
  M[,1] <- M[,ID]
  M[,ID] <- M1
  # the following is the ilr-based transformation for M_1, 
  MT <- matrix(0,n,p-1)
  for (i in 1:n){
    for (j in 1:(p-1)){
      C_1 <- sqrt((p-j)/(p-j+1))
      C_2 <- prod(M[i,(j+1):p]^(1/(p-j)))
      MT[i,j] <- C_1*log(M[i,j]/C_2)
    }
  }
  return(MT)
}