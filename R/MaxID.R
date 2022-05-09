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
#'

#Pr_M = function(P.min, M,
#                ref.geno = NULL,
#                ref.bed = NULL,
#                block.map = NULL,
#                var.name = NULL,
#                output.dir = "./temp/",
#                LD.thresh,
#                K = 1000)
Pr_M = function(M){

  F_min_M = mean(P.min <= M)
  I.m = which(P.min <= M)
  I.m.size = length(I.m)
  P.min.target =  P.min[I.m]
  if(is.null(var.name)){
    var.name.target = I.m
  }else{
    var.name.target =  var.name[I.m]
  }
  names(P.min.target) = var.name.target
  pr.dense = rep(pr.dense.ind(M, F_min_M), I.m.size)

  if(is.null(ref.geno) & is.null(ref.bed) ){
    return(prod(pr.dense))
  }
  #------------------ read the genotype file ------------------------
  if(!is.null(ref.bed)){
    ref.geno = snpStats::read.plink( bed = ref.bed,
                                     select.snps = var.name.target )
    ref.geno = as(ref.geno$genotypes, "numeric")
  }
  #--------------------- define the block list -------------------
  if(is.null(block.map)){
    block.list = Find.plink.block(snplist = var.name.target,
                                  ref.bed = ref.bed,
                                  output.dir = output.dir)
  }else{
    block.list = split(var.name.target, block.map[I.m])
  }
  n.block = length(block.list)
  n.eff = n.eff.vec(ref.geno, block.list)
  return(( prod(pr.dense) )^(n.eff/I.m.size))
}


MinID = function(P1, P2,type ="FWER", alpha = 0.05,
                 ref.geno = NULL,
                 ref.bed = NULL,
                 block.map = NULL,
                 var.name = NULL,
                 output.dir = "./temp/",
                 LD.thresh,
                 K = 1000){

    #-------------------------------------------
    # If the criteria is to control the FWER
    #--------------------------------------------
    P.min = pmin(P1, P2)
    P.max = pmax(P1, P2)
    M = min(P.max)
    M.choice = sort(unique(c(P.max) ), decreasing = FALSE)[1:500]

    if(type == "FWER"){  # minimum p value

      if((1 - Pr_M(M.choice[1])) > alpha ){return(list(S = NULL, tau = NULL))}else{

        pr.m.choice = (sapply(M.choice, Pr_M_prod, P.min = P.min))^K

        tau = max( M.choice[ (1 - pr.m.choice) <= alpha])

        S =  which(P.max <= tau)

      return(list(S = S, tau = tau))
      }
    }

    #------------------------------
    # If the criteria is FDR
    #------------------------------
    if(type == "FDR"){
      p.choice = sort(unique(c(P.max, P.min) ),
                      decreasing = FALSE)[1:1000]
      G = function(t) {sum( P.max < t )/length(P.max) }
      G <- Vectorize(G, SIMPLIFY = T)
      FDR.t = CF(p.choice)/max(1/p, G(p.choice) )
      if(sum( FDR.t <= alpha ) > 0){
        tau = max(p.choice[FDR.t <= alpha])
        S = which(P.max <= tau)
        return(list(S = S, tau = tau))
      }else{
        return(list(S = 0, tau = 0))
        }
    }
}



Pr_M_prod = function( M, P.min){
  F_min_M = mean(P.min <= M)
  I.m = which(P.min <= M)
  I.m.size = length(I.m)
  P.min.target =  P.min[I.m]
  pr.dense = rep(pr.dense.ind(M, F_min_M), I.m.size)
  return(prod(pr.dense))
}


Pv_M_prod = function( M, P.min, Kappa){
  Pv = 1 - Pr_M_prod( M, P.min)^(Kappa)
  return(Pv)
}

MinID.simple = function(P1, P2, Kappa, type = "FWER", alpha = 0.05){
  P.min = pmin(P1, P2)
  P.max = pmax(P1, P2)
  M = min(P.max)
  M.choice = sort(unique(c(P.max) ), decreasing = FALSE)[1:500]
  if(type == "FWER"){  # minimum p value

    if( Pv_M_prod(M.choice[1], P.min, Kappa)  > alpha ){return(list(S = NULL, tau = NULL))}else{

      pv = sapply(M.choice, Pv_M_prod, P.min = P.min, Kappa = Kappa)

      tau = max( M.choice[pv <= alpha])

      S =  which(P.max <= tau)

      return(list(S = S, tau = tau))
    }
  }
}


