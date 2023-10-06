#' Detect function for the genotype data
#'
#' @param P.max
#' @param P.min
#' @param M
#'
#' @return detection pvalue
#'
#' @export
#'
#' @examples

MinMaxP.detect = function(P.max, P.min, M){
  Tstat = M
  F_min_M = mean(P.min <= M)
  I.m = which(P.min <= M)
  I.m.size = length(I.m)
  P.min.target =  P.min[I.m]
  C_0_m = Cmax_null(t = M, P.min = P.min)
  Pv = 1 - exp(- I.m.size * (M * ( F_min_M - M  ))/ ( (1 - M) * F_min_M ) )
  PV2 = 1 - exp(- I.m.size * C_0_m )
  return(c(Pv, PV2) )
}

#fdr_value = function(x){ Fmax_null(t = x, P.min = P.min)/(1/d+mean(P.max < x))  }
#fwer_value = function(x){ Cmax_null(t = x, P.min = P.min) * sum(P.min < x)  }

#' Detect function for the genotype data
#'
#' @param P.max
#' @param P.min
#' @param alpha # significance value
#' @param method # FWER and FDR
#'
#' @return detection pvalue
#'
#' @export
#'
#' @examples

MinMaxP.discov= function(P.max, P.min, method = c("FWER", "FDR"), alpha = 0.05){
  #--------- decide the range of the cut-off  ---------
  # step 1: grid search to determine the initial interval to aovid unnecessary computatation
  d = length(P.min)
  M = min(P.max)
  M.choice.0 = sort(exp( log(M) * c(1:10) * 0.1 ))
  C_0_m = Cmax_null(t = M, P.min = P.min)
  if(method == "FDR"){ # fdr_value
    citeria_fun = function(x){ Fmax_null(t = x, P.min = P.min)/(1/d+mean(P.max < x))  }
  }
  if(method == "FWER"){ #fwer_value
    citeria_fun = function(x){ Cmax_null(t = x, P.min = P.min) * sum(P.min < x)  }
  }
  citeria_value = sapply(M.choice.0, citeria_fun)
  #fdr_init = sapply(M.choice.0, function(x){ Fmax_null(t = x, P.min = P.min)/(1/d+mean(P.max < x))  } )
  #fwer_init = sapply(M.choice.0, function(x){ Cmax_null(t = x, P.min = P.min) * sum(P.min < x)  } )

  # step 2:choose Pmax based on the fdr_init
  k = min(which(citeria_value > alpha))
  M.lower = M.choice.0[k]
  M.upper = max(P.max[P.max <= M.choice.0[k-1]] )
  M.choice = unique(P.max[P.max < M.lower])
  M.choice = sort(M.choice[M.choice >= M.upper ]) # now we had a smaller set of M.choice

  # step 3: elaborate the choice
  pv = sapply(M.choice, citeria_fun)
  tau = max( M.choice[pv <= alpha])
  S = which(P.max <= tau)
  return(list(S = S, tau = tau))
}

# utilities function
Cmax_null = function(t, P.min){
  Fmin_t = mean(P.min <= t)
  return( t * ( Fmin_t -  t ) / ( (1 - t)* Fmin_t )  )
}

Fmax_null = function(t, P.min){
  Fmin_t = mean(P.min <= t)
  return( t * ( Fmin_t -  t ) / ( (1 - t))  )
}




#
# MinMaxP.geno = function(P1, P2,
#                    ref.geno = NULL,
#                    ref.bed = NULL,
#                    block.map = NULL,
#                    var.name = NULL,
#                    output.dir = "./temp/",
#                    block.thresh = 0.995,
#                    K = 1000,
#                    exact = F){
#   # the input is p-value sequence
#   # geno is a bigsnpr object without missing values
#
#   P.min = pmin(P1, P2)
#   P.max = pmax(P1, P2)
#   M = min(P.max)
#   Tstat = M
#   F_min_M = mean(P.min <= M)
#   I.m = which(P.min <= M)
#   I.m.size = length(I.m)
#   P.min.target =  P.min[I.m]
#   if(is.null(var.name)){
#     var.name.target = I.m
#   }else{
#     var.name.target =  var.name[I.m]
#   }
#   names(P.min.target) = var.name.target
#   pr.sparse = pr.sparse.ind(M, P.min.target)
#   #pr.dense = rep(pr.dense.ind(M, F_min_M), I.m.size)
#   log.pr.dense = log(1 - M)/F_min_M + ((1- F_min_M)/F_min_M) * (M - M^2/2)/(1 - M)^2
#
#   if(is.null(ref.geno) & is.null(ref.bed) ){
#     cat("Ignoring the dependence")
#     return(list(M = M, pv.sparse = 1 - prod(pr.sparse),
#                 pv.dense = 1 - exp(I.m.size * log.pr.dense) ))
#   }
#   #------------------ read the genotype file ------------------------
#   if(!is.null(ref.geno)){
#     loc.target = match(var.name.target, ref.geno$map$marker.ID)
#     ref.geno = ref.geno$genotypes[,loc.target] %>% as.matrix
#     colnames(ref.geno) = var.name.target
#   }else{
#
#   if(!is.null(ref.bed)){
#      ref.geno = snpStats::read.plink( bed = ref.bed,
#                                       select.snps = var.name.target )
#      ref.geno = as(ref.geno$genotypes, "numeric")
#     }
#   }
#   #--------------------- define the block list -------------------
#   if(is.null(block.map)){
#   block.list = Find.plink.block(snplist = var.name.target,
#                                 ref.bed = ref.bed,
#                                 output.dir = output.dir)
#   }else{
#     block.list = split(var.name.target, block.map[I.m])
#   }
#   n.block = length(block.list)
#   n.eff = n.eff.vec(ref.geno, block.list, block.thresh = block.thresh)
#   p.value = 1 - exp(n.eff * log.pr.dense)
#   return(list(M = Tstat, p = p.value, pr.dense = exp(I.m.size * log.pr.dense), I.m.size = (I.m.size),
#               n.block = n.block, n.eff = n.eff))
# }
#
# Adj.func = function(P1, P2){
#   P.min = pmin(P1, P2)
#   P.max = pmax(P1, P2)
#   M = min(P.max)
#   I.m = which(P.min <= M)
#   P.min.target =  P.min[I.m]
#   #I.1 = which(mu != 0)
#   #I.0 = which(mu == 0)
#
#   pr.approx = (1-M)/(1-P.min.target)
#   #F_min = pp.to.z(mu = mu, px =M) + M - M*(pp.to.z(mu = mu, px =M))
#   F_min.est = rep(mean(P.min < M), length(P.min.target))
#   #A.1 = sum((1 - F_min) /F_min)
#   A.1 = sum((1 - F_min.est) /F_min.est)
#   A.2 = -log(1 - M) - (M - M^2/2)/(1 - M)^2
#   A.3 = sum(log(1 - P.min.target))
#
#   A = A.1*A.2 - A.3
#
#   return( exp(-A) )
#
# }
#


