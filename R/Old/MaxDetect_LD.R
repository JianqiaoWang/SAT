#' Detect function
#'
#' @param P1
#' @param P2
#' @param K
#' @param ref.geno
#'
#' @return
#' @export
#'
#' @examples
MinDetect_LD = function(P1, P2, K = 1000, ref.geno = NULL){
  # the input is p-value sequence
  vec = vector()
  P.min = pmin(P1, P2)
  P.max = pmax(P1, P2)
  Tstat = min(P.max)
  target.loci = which(P.min < Tstat)
  P.min.target =  P.min[target.loci]
  ref.geno = scale(ref.geno[,target.loci])
  vec= vector()
  for(i in 1:K){
    Y = rnorm(nrow(ref.geno), 0, 1); # simulate the probability from the null
    Z <- crossprod(Y, ref.geno)/sqrt(nrow(ref.geno));
    p.z =  2* (1 - pnorm( abs(Z)))
    p.tilde = (1 - P.min.target) * p.z + P.min.target;
    vec = c(vec, min(p.tilde));
  }

  p.value = sum(vec <= Tstat)/length(vec)
  return(list(M = Tstat, p = p.value))
}

