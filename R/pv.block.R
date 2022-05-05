pr.sparse.ind = function(M, P.min){
  # pr( P.max > Tstat |P.min) for sparse effects
  pr = (1 - M)/(1 - P.min);
  return(pr)
}

pr.dense.ind = function(M, F_min_M){
  # pr( P.max > Tstat |P.min) for dense effects
  pr = exp(log(1 - M)/F_min_M + ((1- F_min_M)/F_min_M) * (M - M^2/2)/(1 - M)^2 )
  return(pr)
}

#-------------------------------------------------------------------
# pr( P.max > Tstat |P.min) for a given block
# Monte-Carlo simulation
# Specify the covariance matrix and
#--------------------------------------------------------------------
pr.block.sparse.dep = function(Tstat, P.min, ref.geno, exact = F){
  block.size = length(P.min)
  if(block.size == 1){
    pr = (1 - Tstat)/(1 - P.min); # pr( P.max > Tstat |P.min)
  }else{
    # for the block with more than one variants, we consider the monte-carlo method to simulate
    # min P_max within block i; after transforming it to the z score.

    if(exact == T){
      cond.pv = (Tstat -P.min)/(1 - P.min);
      abs.z = qnorm(cond.pv/2, lower.tail = F)
      mean <- rep(0, block.size)
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

pv.block.add.groupinfo = function(Tstat, P.min, ref.geno, which.min.index, exact = F){
  # pr( P.max > Tstat |P.min) for a given block
  if(length(P.min) == 1){
    pr = (1 - Tstat)/(1 - P.min); # pr( P.max > Tstat |P.min)
  }else{
    if(exact == T){
      cond.pv = (Tstat -P.min)/(1 - P.min);
      abs.z = qnorm(cond.pv/2, lower.tail = F)
      n <- length(P.min)
      mean <- rep(0, n)
      corr <- cor(ref.geno)

      if(length(unique(which.min.index) ) == 1){
        pr <- pmvnorm(-abs.z, abs.z, mean, corr, algorithm = GenzBretz(abseps = 0.000001))
        pr = pr[1]
      }

      if(length(unique(which.min.index) ) == 2){
      group1 = which(which.min.index == 1)
      group2 = which(which.min.index == 2)
      pr.seq1 <- pmvnorm(-abs.z[group1], abs.z[group1],
                         mean[group1], sigma = corr[group1, group1], algorithm = GenzBretz(abseps = 0.000001))
      pr.seq2 <- pmvnorm(-abs.z[group2], abs.z[group2],
                         mean[group2], sigma = corr[group2, group2], algorithm = GenzBretz(abseps = 0.000001))
      pr = pr.seq1[1] * pr.seq2[1]
      }

    }else{
      # for the block with more than one variants, we consider the monte-carlo method to simulate
      # min P_max with in block i;
      vec = t(replicate(1000,simnull( P.min = P.min,  ref.geno = ref.geno)))
      # then calculate the probability of P(min P_max >  Tstat |P.min)
      pr = sum(vec >= Tstat)/length(vec)
    }

  }
  return(pr)

}

pv.block2.add.groupinfo = function(Tstat, P.min, ref.geno, which.min.index, exact = F){
  # pr( P.max > Tstat |P.min) for a given block
  if(length(P.min) == 1){
    pr = (1 - Tstat)/(1 - P.min); # pr( P.max > Tstat |P.min)
  }else{
    # for the block with more than one variants, we consider the monte-carlo method to simulate
    # min P_max with in block i;
    if(exact == T){
      #conditional probability of P.max larger than P.min
      abs.z.min = qnorm(P.min/2, lower.tail = F)
      n <- length(P.min)
      mean <- rep(0, n)
      corr <- cor(ref.geno)
      abs.z.Tstat = rep(qnorm(Tstat/2, lower.tail = F),n)

      cond.pr = function(abs.z.Tstat, abs.z.min, mean, corr){
        pr.larger.p.min <- pmvnorm(lower = -abs.z.min, upper = abs.z.min,
                                   mean =  mean, sigma = corr, algorithm = GenzBretz(abseps = 0.000001))
        pr.larger.p.Tstat <- pmvnorm(lower =-abs.z.Tstat, upper = abs.z.Tstat, mean = mean,
                                    sigma = corr,algorithm = GenzBretz(abseps = 0.000001))
        pr = pr.larger.p.Tstat[1]/pr.larger.p.min[1]
        return(pr)
      }


      if(length(unique(which.min.index) ) == 1){
        pr = cond.pr(abs.z.Tstat, abs.z.min, mean, corr)
      }

      if(length(unique(which.min.index) ) == 2){
        group1 = which(which.min.index == 1) # symmetric for group 1 or group 2
        group2 = which(which.min.index == 2)
        pr.seq1 <- cond.pr(abs.z.Tstat = abs.z.Tstat[group1], abs.z.min =  abs.z.min[group1],
                           mean = mean[group1], corr = corr[group1, group1])
        pr.seq2 <- cond.pr(abs.z.Tstat =  abs.z.Tstat[group2], abs.z.min = abs.z.min[group2],
                           mean = mean[group2], corr = corr[group2, group2])
        pr = pr.seq1 * pr.seq2
      }



    }else{
    vec = t(replicate(1000,simnull( P.min = P.min,  ref.geno = ref.geno)))
    # then calculate the probability of P(min P_max >  Tstat |P.min)
    pr = sum(vec >= Tstat)/length(vec)
    }

  }
  return(pr)

}


pv.block2 = function(Tstat, P.min, ref.geno,  exact = F){
  # pr( P.max > Tstat |P.min) for a given block; consider the location of the p-values
  if(length(P.min) == 1){
    pr = (1 - Tstat)/(1 - P.min); # pr( P.max > Tstat |P.min)
  }else{
    # for the block with more than one variants, we consider the monte-carlo method to simulate
    # min P_max with in block i;
    if(exact == T){
      #conditional probability of P.max larger than P.min
      abs.z.min = qnorm(P.min/2, lower.tail = F)
      n <- length(P.min)
      mean <- rep(0, n)
      corr <- cor(ref.geno)
      pr.larger.p.min <- pmvnorm(-abs.z.min, abs.z.min,
                                 mean, corr, algorithm = GenzBretz(abseps = 0.000001))
      abs.z.Tstat = rep(qnorm(Tstat/2, lower.tail = F),n)
      pr.larger.p.Tstat <- pmvnorm(-abs.z.Tstat, abs.z.Tstat, mean, corr,algorithm = GenzBretz(abseps = 0.000001))
      pr = pr.larger.p.Tstat/pr.larger.p.min
    }else{
      vec = t(replicate(1000,simnull( P.min = P.min,  ref.geno = ref.geno)))
      # then calculate the probability of P(min P_max >  Tstat |P.min)
      pr = sum(vec >= Tstat)/length(vec)
    }

  }
  return(pr)

}


pv.block3.add.groupinfo = function(Tstat, P.min, ref.geno, which.min.index, exact = F){
  # pr( P.max > Tstat |P.min) for a given block
  if(length(P.min) == 1){
    pr = (1 - Tstat)/(1 - P.min); # pr( P.max > Tstat |P.min)
  }else{
    # for the block with more than one variants, we consider the monte-carlo method to simulate
    # min P_max with in block i;
    if(exact == T){
      #conditional probability of P.max larger than P.min
      abs.z.min = qnorm(P.min/2, lower.tail = F)
      n <- length(P.min)
      mean <- rep(0, n)
      corr <- cor(ref.geno)
      abs.z.Tstat = rep(qnorm(Tstat/2, lower.tail = F),n)

      cond.pr = function(abs.z.Tstat, abs.z.min, mean, corr){
        pr.larger.p.min <- pmvnorm(lower = -abs.z.min, upper = abs.z.min,
                                   mean =  mean, sigma = corr, algorithm = GenzBretz(abseps = 0.000001))
        pr.larger.p.Tstat <- pmvnorm(lower =-abs.z.Tstat, upper = abs.z.Tstat, mean = mean,
                                     sigma = corr,algorithm = GenzBretz(abseps = 0.000001))
        pr = pr.larger.p.Tstat[1]/pr.larger.p.min[1]
        return(pr)
      }


      if(length(unique(which.min.index) ) == 1){
        pr = cond.pr(abs.z.Tstat, abs.z.min, mean, corr)
      }

      if(length(unique(which.min.index) ) == 2){
        group1 = which(which.min.index == 1)
        group2 = which(which.min.index == 2)
        # calculate the conditional distribution for seq 1: that is conditional on the pmin in group 1 ,
        # calculate the distribution of elements of index of group 2;
        conditional.para =
          condMVN(mean = mean, sigma = corr, dependent.ind = group2,
                  given.ind = group1 , X.given = abs.z.min[group1], check.sigma=TRUE)
        conditional.mean = conditional.para$condMean
        conditional.corr = conditional.para$condVar
        pr.seq1 <- cond.pr(abs.z.Tstat = abs.z.Tstat[group2],
                           abs.z.min =  abs.z.min[group2],
                           mean = conditional.mean,
                           corr = conditional.corr)

        conditional.para2 =
          condMVN(mean = mean, sigma = corr, dependent.ind = group1,
                  given.ind = group2 , X.given = abs.z.min[group2], check.sigma=TRUE)
        conditional.mean2 = conditional.para2$condMean
        conditional.corr2 = conditional.para2$condVar

        pr.seq2 <- cond.pr(abs.z.Tstat =  abs.z.Tstat[group1],
                           abs.z.min = abs.z.min[group1],
                           mean = conditional.mean2,
                           corr = conditional.corr2)
        pr = pr.seq1 * pr.seq2
      }



    }else{
      vec = t(replicate(1000,simnull( P.min = P.min,  ref.geno = ref.geno)))
      # then calculate the probability of P(min P_max >  Tstat |P.min)
      pr = sum(vec >= Tstat)/length(vec)
    }

  }
  return(pr)

}



pv.block3 = function(Tstat, P.min, ref.geno){
  # pr( P.max > Tstat |P.min) for a given block
  if(length(P.min) == 1){
    pr = (1 - Tstat)/(1 - P.min); # pr( P.max > Tstat |P.min)
  }else{
    # for the block with more than one variants, we consider the monte-carlo method to simulate
    # min P_max with in block i;
    vec = t(replicate(1000,simnull.3( P.min = P.min,  ref.geno = ref.geno )))
    # then calculate the probability of P(min P_max >  Tstat |P.min)
    pr = sum(vec >= Tstat)/length(vec)

  }
  return(pr)

}




simnull2 = function( P.min, ref.geno){
  # method 1 for simulate the mini-max pvalue in a given block : require the ref.geno is scaled
  Y = rnorm(nrow(ref.geno), 0, 1); # simulate the probability from the null
  r <- cor(Y, ref.geno);
  Z <-  r.to.z(r= r, n = nrow(ref.geno));
  p.z =  2* (1 - pnorm( abs(Z)))
  p.tilde = (1 - P.min) * p.z + P.min;
  return(min(p.tilde))
}


