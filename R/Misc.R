simnull.3 = function( P.min, ref.geno){
  # method 1 for simulate the mini-max pvalue in a given block : require the ref.geno is scaled
  Y = rnorm(nrow(ref.geno), 0, 1); # simulate the probability from the null
  Z <- crossprod(Y, ref.geno)/sqrt(nrow(ref.geno));
  p.z =  2* (1 - pnorm( abs(Z)))
  p.tilde = (1 - P.min) * p.z + P.min;
  return(min(p.tilde))
}


simnull = function( P.min, ref.geno){
  # method 1 for simulate the mini-max pvalue in a given block : require the ref.geno is scaled
  Y = rnorm(nrow(ref.geno), 0, 1); # simulate the probability from the null
  r <- cor(Y, ref.geno);
  Z <-  r.to.z(r= r, n = nrow(ref.geno));
  p.z =  2* (1 - pnorm( abs(Z)))
  p.tilde = (1 - P.min) * p.z + P.min;
  return(min(p.tilde))
}

condMVN = function (mean, sigma, dependent.ind, given.ind, X.given, check.sigma = TRUE)
{
  if (missing(dependent.ind))
    return("You must specify the indices of dependent random variables in `dependent.ind'")
  if (missing(given.ind) & missing(X.given))
    return(list(condMean = mean[dependent.ind], condVar = as.matrix(sigma[dependent.ind,
                                                                          dependent.ind])))
  if (length(given.ind) == 0)
    return(list(condMean = mean[dependent.ind], condVar = as.matrix(sigma[dependent.ind,
                                                                          dependent.ind])))
  if (length(X.given) != length(given.ind))
    stop("lengths of `X.given' and `given.ind' must be same")
  if (check.sigma) {
    if (!isSymmetric(sigma))
      stop("sigma is not a symmetric matrix")
    eigenvalues <- eigen(sigma, only.values = TRUE)$values
    if (any(eigenvalues < 1e-08))
      stop("sigma is not positive-definite")
  }
  B <- sigma[dependent.ind, dependent.ind]
  C <- sigma[dependent.ind, given.ind, drop = FALSE]
  D <- sigma[given.ind, given.ind]
  CDinv <- C %*% solve(D)
  cMu <- c(mean[dependent.ind] + CDinv %*% (X.given - mean[given.ind]))
  cVar <- B - CDinv %*% t(C)
  list(condMean = cMu, condVar = cVar)
}


