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

MinID = function(P1, P2,type ="FWER", alpha = 0.05, eff.ratio = 0.25,
                 ref.geno = NULL,
                 ref.bed = NULL,
                 block.map = NULL,
                 var.name = NULL,
                 output.dir = "./temp/",
                 simple = T,
                 K = 1000){

    P.min = pmin(P1, P2)
    P.max = pmax(P1, P2)
    M = min(P.max)

    if(type == "FWER"){
      #-------------------------------------------
      # If the criteria is to control the FWER
      #--------------------------------------------
      if(Pv_M_prod(M, P.min, eff.ratio) > alpha ){return(list(S = NULL, tau = NULL))}else{

        tau = Search.simple(P.min, P.max, eff.ratio, alpha = alpha)

        if(!simple){ # one more round of updating eff ratio
        Pv.tau = Pv_given_M(M = tau, P.min = P.min, ref.geno = ref.geno, ref.bed = ref.bed,
             block.map = block.map, var.name = var.name, output.dir = output.dir )

        print(Pv.tau)

        eff.ratio = Pv.tau$n.eff/Pv.tau$I.m.size  # update the eff.ratio

        print(eff.ratio)

        tau = Search.simple(P.min = P.min, P.max = P.max, eff.ratio = eff.ratio, alpha = alpha)
        }
        S =  which(P.max <= tau)
        }
      return(list(S = S, snp.name = var.name[S], tau = tau, eff.raio = eff.ratio))
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



log.pr.dense.fun = function(M, F_min_M){
  log.pr.dense = log(1 - M)/F_min_M + ((1- F_min_M)/F_min_M) * (M - M^2/2)/(1 - M)^2
  return(log.pr.dense)
}


Pv_M_prod = function( M, P.min, eff.ratio, all.len = length(P.min) ){
  F_min_M = sum(P.min <= M)/all.len
  I.m.size = sum(P.min <= M)
  log.pr.dense = log(1 - M)/F_min_M + ((1- F_min_M)/F_min_M) * (M - M^2/2)/(1 - M)^2
  Pv = 1 - exp( I.m.size * eff.ratio *  log.pr.dense )
  return(Pv)
}


Search.simple = function(P.min, P.max, eff.ratio, type = "FWER", alpha = 0.05){
  # search with fixed eff.ratio
  # return tau, the max cut-off value satisfying the < alpha
  M = min(P.max)
  if(type == "FWER"){
      #--------- decide the the range of M choice ---------
      all.len = length(P.min)
      M.choice.0 = sort(exp( log(M) * c(1:10) * 0.1 ))
      for(k in 2:length(M.choice.0)){
        if(Pv_M_prod(M.choice.0[k], P.min, eff.ratio) > alpha) break
      }

      M.lower = M.choice.0[k]
      M.upper = max(P.max[P.max <= M.choice.0[k-1]] )
      M.choice = unique(P.max[P.max < M.lower])
      M.choice = sort(M.choice[M.choice >= M.upper ])
      P.min = P.min[P.min < max(M.choice) ]
      pv = sapply(M.choice, Pv_M_prod, P.min = P.min, eff.ratio = eff.ratio, all.len = all.len)
      print(pv)
      tau = max( M.choice[pv <= alpha])
      return(tau)
  }
}


Pv_given_M = function(M, P.min,
                ref.geno = NULL,
                ref.bed = NULL,
                block.map = NULL,
                var.name = NULL,
                output.dir = "./temp/",
                block.thresh = 0.995,
                K = 1000,
                exact = F){
  # the input is p-value sequence
  # geno is a bigsnpr object without missing values
  Tstat = M
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
  pr.sparse = pr.sparse.ind(M, P.min.target)
  log.pr.dense = log(1 - M)/F_min_M + ((1- F_min_M)/F_min_M) * (M - M^2/2)/(1 - M)^2

  if(is.null(ref.geno) & is.null(ref.bed) ){
    cat("Ignoring the dependence")
    return(list(M = M, pv.sparse = 1 - prod(pr.sparse),
                pv.dense = 1 - exp(I.m.size * log.pr.dense) ))
  }
  #------------------ read the genotype file ------------------------
  if(!is.null(ref.geno)){
    loc.target = match(var.name.target, ref.geno$map$marker.ID)
    ref.geno = ref.geno$genotypes[,loc.target] %>% as.matrix
    colnames(ref.geno) = var.name.target
  }else{

    if(!is.null(ref.bed)){
      ref.geno = snpStats::read.plink( bed = ref.bed,
                                       select.snps = var.name.target )
      ref.geno = as(ref.geno$genotypes, "numeric")
    }
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
  n.eff = n.eff.vec(ref.geno, block.list, block.thresh = block.thresh)[2]
  p.value = 1 - exp(n.eff * log.pr.dense)
  return(list(M = Tstat, p = p.value, pr.dense = exp(I.m.size * log.pr.dense), I.m.size = (I.m.size),
              n.block = n.block, n.eff = n.eff))
}

