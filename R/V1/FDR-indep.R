#' Identification of sharing association
#'
#' @param P1 the input is p-value sequence
#' @param P2
#' @param K
#' @param F.null null distribution
#' @param type Criteria for the identification
#' @param alpha error rate level
#'
#' @return return is S and tau, where tau is the p-value threshold
#' @export
#'
#' @examples
#'
#'
#'discov

identification = function(P1, P2,type ="FDR", alpha = 0.05, eff.ratio = 0.25,
                 ref.geno = NULL,
                 ref.bed = NULL,
                 block.map = NULL,
                 var.name = NULL,
                 output.dir = "./temp/",
                 simple = T,
                 K = 1000){
  
  P.min = pmin(P1, P2)
  P.max = pmax(P1, P2)
  Fmin.bar = function(t){ mean(P.min < t) }
  Fmax.bar = function(t){ mean(P.max < t) }
  FDR = function(t){ t* ( Fmin.bar(t) - t )/ ( (1 - t)* Fmax.bar(t) ) }
  #M = min(P.max)
  p.choice = sort(unique(c(P.max, P.min) ),
                  decreasing = FALSE)[1:1000]
  
}

