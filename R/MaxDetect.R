#' Detect function
#'
#' @param P1 
#' @param P2 
#' @param K 
#' @param F.null 
#'
#' @return
#' @export
#'
#' @examples
MinDetect = function(P1, P2, K = 200, F.null = NULL){
  # the input is p-value sequence 
  vec = vector()
  P.min = pmin(P1, P2)
  P.max = pmax(P1, P2)
  n= length(P.min)
  Tstat = min(P.max)
  # conditional on the P.min, the distribution of P.max is given by 
  if(is.null(F.null)){
    pr = (1 - Tstat)/(1 - P.min); # pr( P.max > Tstat |P.min)
    pr[pr>1] = 1; pr[pr < 0] = 0; # truncate it to 0, 1 interval
    p.value = 1 - prod(pr); # calculate the p value for the observed stat
    }else{
    for(i in 1:K){
    p.z = F.null(); # simulate the probability from the null 
    p.tilde = (1 - P.min) * p.z + P.min;
    vec = c(vec, min(p.tilde) );
    }
      p.value = sum(vec <= Tstat)/length(vec)
  }
  return(list(M = Tstat, p = p.value))
}
# 
# MaxDetect.GWAS.Z = function(P1, P2, K = 200, Sigma = NULL){
#   # the input is p-value sequence 
#   vec = vector()
#   P.min = pmin(P1, P2)
#   P.max = pmax(P1, P2)
#   n= length(P.min)
#   Tstat = min(P.max)
#   # conditional on the P.min, the distribution of P.max is given by 
#   if(is.null(Sigma)){
#     pr = (1 - Tstat)/(1 - P.min) # pr( P.max > Tstat |P.min)
#     pr[pr>1] = 1; pr[pr < 0] = 0; # truncate it to 0, 1 interval
#     p.value = 1 - prod(pr) # calculate the p value for the observed stat
#   }
#   return(list(M = Tstat, p = p.value))
# }
# 
# 
