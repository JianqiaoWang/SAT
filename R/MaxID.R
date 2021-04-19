#' Identification of sharing association
#'
#' @param P1
#' @param P2
#' @param K
#' @param F.null null distribution
#' @param type Criteria for the identification
#' @param alpha error rate level
#'
#' @return
#' @export
#'
#' @examples
MinID = function(P1, P2, K = 200, F.null = NULL, type, alpha){
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
    };CF <- Vectorize(CF, SIMPLIFY = T)

    if(type == "FWER"){  # minimum p value
      t.choice = sort(unique(P.max), decreasing = FALSE)[1:1000]
      if(CF(Tstat)){return(list(S = 0, tau = 0))}else{
        tau = min(t.choice[CF(t.choice) <= alpha])
        S =  which(P.max <= tau)
      return(list(S = S, tau = tau))
      }
    }

    if(type == "FDR"){
      t.choice = sort(unique(c(P.max, P.min) ),
                      decreasing = FALSE)[1:1000]
      G = function(t) {sum( P.max < t )/length(P.max) }
      G <- Vectorize(G, SIMPLIFY = T)
      FDR.t = CF(t.choice)/max(1/d, G(t.choice) )
      if(sum( FDR.t <= alpha ) > 0){
        tau = min(t.choice[FDR.t <= alpha])
        S = which(P.max >= t_opt)
        return(list(S = S, tau = tau))
      }else{
        return(list(S = 0, tau = 0))
        }
    }
}

