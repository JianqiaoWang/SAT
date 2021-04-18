library(mvtnorm)
library(Matrix)
bdiag_m <- function(lmat) {
  ## Copyright (C) 2016 Martin Maechler, ETH Zurich
  if(!length(lmat)) return(new("dgCMatrix"))
  stopifnot(is.list(lmat), is.matrix(lmat[[1]]),
            (k <- (d <- dim(lmat[[1]]))[1]) == d[2], # k x k
            all(vapply(lmat, dim, integer(2)) == k)) # all of them
  N <- length(lmat)
  if(N * k > .Machine$integer.max)
    stop("resulting matrix too large; would be  M x M, with M=", N*k)
  M <- as.integer(N * k)
  ## result: an   M x M  matrix
  new("dgCMatrix", Dim = c(M,M),
      ## 'i :' maybe there's a faster way (w/o matrix indexing), but elegant?
      i = as.vector(matrix(0L:(M-1L), nrow=k)[, rep(seq_len(N), each=k)]),
      p = k * 0L:M,
      x = as.double(unlist(lmat, recursive=FALSE, use.names=FALSE)))
}

#------------------------------------- option 1
GenZ1 = function(SIGMA){
Z = lapply(1:nblock, function(i) {
  mvtnorm::rmvnorm(1,  sigma = SIGMA[[i]] )
}) %>% unlist()
}

#------------------------------------- option 2
GenZ2 = function(SIGMA.SQRT){
Z = SIGMA.SQRT %*% rnorm(dim(SIGMA.SQRT)[1] ) %>% as.vector()
}
#------------------------------------ Distribution function of |Z| under the null
# ----------------------------------- the ouput is p value 
F.null.Z = function(SIGMA.SQRT){
  Z = SIGMA.SQRT %*% rnorm(dim(SIGMA.SQRT)[1] ) %>% as.vector()
  p.z = 2*(1 - pnorm( abs(z)))
}


#microbenchmark::microbenchmark(a = GenZ2(SIGMA.SQRT), b = GenZ1(SIGMA), times = 5)
