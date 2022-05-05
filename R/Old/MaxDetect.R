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
    vec = c(vec, min(p.tilde));
    }
      p.value = sum(vec <= Tstat)/length(vec)
  }
  return(list(M = Tstat, p = p.value))
}

Adj.func = function(P1, P2){
  P.min = pmin(P1, P2)
  P.max = pmax(P1, P2)
  M = min(P.max)
  target.loci = which(P.min <= M)
  P.min.target =  P.min[target.loci]
  #I.1 = which(mu != 0)
  #I.0 = which(mu == 0)

  pr.approx = (1-M)/(1-P.min.target)
  #F_min = pp.to.z(mu = mu, px =M) + M - M*(pp.to.z(mu = mu, px =M))
  F_min.est = rep(mean(P.min < M), length(P.min.target))
  #A.1 = sum((1 - F_min) /F_min)
  A.1 = sum((1 - F_min.est) /F_min.est)
  A.2 = -log(1 - M) - (M - M^2/2)/(1 - M)^2
  A.3 = sum(log(1 - P.min.target))

  A = A.1*A.2 - A.3

  return( exp(-A) )

}


MinDetect.adj = function(P1, P2){
  P.min = pmin(P1, P2)
  P.max = pmax(P1, P2)
  M = min(P.max)
  target.loci = which(P.min <= M)
  P.min.target =  P.min[target.loci]
  pr.approx = (1-M)/(1-P.min.target)
  Adj = Adj.func(P1, P2)
  p.value = 1 - prod(pr.approx)* Adj
  return(list(M = M, p = p.value))

}

MinDetect_Dep = function(P1, P2, K = 1000, corr = NULL, block.map = NULL, LD.thresh, exact = F){
  # the input is p-value sequence
  # When LD.thresh is 1, it degrades to the independence case
  # When LD.thresh is 0, it degrades to the total dependence case
  vec = vector()
  P.min = pmin(P1, P2)
  P.max = pmax(P1, P2)
  Tstat = min(P.max)
  target.loci = which(P.min <= Tstat)
  P.min.target =  P.min[target.loci]
  corr = corr[target.loci,target.loci] %>% as.matrix()

  #--------------------- define the block list -------------------
  if(is.null(block.map) ){
    adj.cor.mat = sapply(1:(ncol(ref.geno)-1), function(x){
      cor(ref.geno[,x], ref.geno[,x+1], use = "pairwise.complete.obs")
    } )
    block.list =  GenerateBlock(adj.cor.mat, rho = LD.thresh)
  }else{

    block.list =  split(1:length(target.loci), block.map[target.loci])

  }
  #-----------------------------------------------------------------
  #block.list.05 =  GenerateBlock(adj.cor.mat, rho = LD.thresh)
  #block.list.07 =  GenerateBlock(adj.cor.mat, rho = LD.thresh)

  pr = sapply(block.list, function(x){ # generate block list  probability
    if(length(x) == 1){
      pr = (1 - Tstat)/(1 - P.min.target[x]); # pr( P.max > Tstat |P.min)
    }else{
      pr = pv.block.basic(Tstat,
                    P.min = P.min.target[x],
                    corr = corr[x,x], exact = T) # pr( P.max > Tstat |P.min) within a block
    }
    return(pr)
  })  # pr( P.max > Tstat |P.min)

  p.value = 1 - prod(pr) ; # calculate the p value for the observed stat

  return(list(M = Tstat, p = p.value))
}

MinDetect_Dep.adj = function(P1, P2, K = 1000, corr = NULL, block.map = NULL, LD.thresh, exact = F){
  # the input is p-value sequence
  # When LD.thresh is 1, it degrades to the independence case
  # When LD.thresh is 0, it degrades to the total dependence case
  vec = vector()
  P.min = pmin(P1, P2)
  P.max = pmax(P1, P2)
  Tstat = min(P.max)
  target.loci = which(P.min <= Tstat)
  P.min.target =  P.min[target.loci]
  corr = corr[target.loci,target.loci] %>% as.matrix()
  M = Tstat
  #--------------------- define the block list -------------------
  if(is.null(block.map) ){
    adj.cor.mat = sapply(1:(ncol(ref.geno)-1), function(x){
      cor(ref.geno[,x], ref.geno[,x+1], use = "pairwise.complete.obs")
    } )
    block.list =  GenerateBlock(adj.cor.mat, rho = LD.thresh)
  }else{

    block.list =  split(1:length(target.loci), block.map[target.loci])

  }
  #-----------------------------------------------------------------
  pr = sapply(block.list, function(x){ # generate block list  probability
    if(length(x) == 1){
      pr = (1 - Tstat)/(1 - P.min.target[x]); # pr( P.max > Tstat |P.min)
    }else{
      pr = pv.block.basic(Tstat,
                          P.min = P.min.target[x],
                          corr = corr[x,x], exact = T) # pr( P.max > Tstat |P.min) within a block
    }
    return(pr)
  })  # pr( P.max > Tstat |P.min)
 #------------------------------------------------------------------
  p.value = 1 - prod(pr) *Adj.func(P1, P2); # calculate the p value for the observed stat

  return(list(M = Tstat, p = p.value))
}

pv.block.basic = function(Tstat, P.min, corr, exact = F){
  if(length(P.min) == 1){
    pr = (1 - Tstat)/(1 - P.min); # pr( P.max > Tstat |P.min)
  }else{
    # for the block with more than one variants, we consider the monte-carlo method to simulate
    # min P_max with in block i;
    if(exact == T){
      cond.pv = (Tstat -P.min)/(1 - P.min);
      abs.z = qnorm(cond.pv/2, lower.tail = F)
      n <- length(P.min)
      mean <- rep(0, n)
      pr <- pmvnorm(-abs.z, abs.z, mean, corr, algorithm = GenzBretz(abseps = 0.000001))
      pr = pr[1]
    }else{
      vec = t(replicate(1000,simnull( P.min = P.min,  ref.geno = ref.geno)))
      # then calculate the probability of P(min P_max >  Tstat |P.min)
      pr = sum(vec >= Tstat)/length(vec)
    }

  }
  return(pr)

}

