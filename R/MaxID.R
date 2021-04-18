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
MinDetect = function(P1, P2, K = 200, F.null = NULL, type, alpha){
  # the input is p-value sequence 
  vec = vector()
  P.min = pmin(P1, P2)
  P.max = pmax(P1, P2)
  n= length(P.min)
  Tstat = min(P.max)
  # conditional on the P.min, the distribution of P.max is given by 
    CF = function(t){
    if(is.null(F.null)){
    pr = (1 - t)/(1 - P.min); # pr( P.max > Tstat |P.min)
    pr[pr>1] = 1; pr[pr < 0] = 0; # truncate it to 0, 1 interval
    p.value = 1 - prod(pr); # calculate the p value for the observed stat
    }else{
    return(1 - F.null(t))
    }
    }
    if(type == "FWER"){  # minimum p value
      t.choice = sort(unique(P.max), decreasing = TRUE)[1:1000]
      tau = min(t.choice[CF(t.choice) <= alpha])
      S =  which(P.min <= tau)
      return(list(S = S, tau = tau))
    }
    if(type == "FDR"){
      t.choice = sort(unique(c(P.max, P.min) ), decreasing = TRUE)[1:1000]
      G = function(t) {sum( T.min > t )/length(T.min) }
      FDR.t = CF(t.choice)/max(1/d, G(t.choice) )
      if(sum( FDR.t <= alpha ) > 0){
        tau = min(t.choice[FDR.t <= alpha])
        S = which(T.min >= t_opt)
      }else{
        return(0)
      }
    }
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
