pv.block = function(Tstat, P.min, ref.geno, exact = F){
  # pr( P.max > Tstat |P.min) for a given block
  if(length(P.min) == 1){
    pr = (1 - Tstat)/(1 - P.min); # pr( P.max > Tstat |P.min)
  }else{
    # for the block with more than one variants, we consider the different method to simulate
    # min P_max within the block i;
    if(exact == T){
      cond.pv = (Tstat -P.min)/(1 - P.min);
      abs.z = qnorm(cond.pv/2, lower.tail = F)
      n <- length(P.min)
      mean <- rep(0, n)
      corr <- cor(ref.geno)
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

p.func = function(P1, P2){
  P.min = pmin(P1, P2)
  P.max = pmax(P1, P2)
  Tstat = min(P.max)
  M = Tstat
  target.loci = which(P.min <= M)
  P.min.target =  P.min[target.loci]
  target.loci = which(P.min <= Tstat)
  F_min_M = mean(P.min <= M)
  pr.block =  rep(exp(log(1 - M)/F_min_M + ((1- F_min_M)/F_min_M) * (M - M^2)/(1 - M)^2 ),
                  length(target.loci))

}


