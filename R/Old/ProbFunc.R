Prob.func = function(P1, P2){

P.min = pmin(P1, P2)
P.max = pmax(P1, P2)
Tstat = min(P.max)
target.loci = which(P.min <= Tstat)
P.min.target =  P.min[target.loci]

# pr(p0 > M | p0 > p_min)
PMAT = cbind(P1[target.loci], P2[target.loci])
group.which.min.index = apply(PMAT, 1, which.min )
cond.pv = (1 - Tstat)/(1 - P.min.target)
pr.p0.geq.Tstat.given.pmin.leq.Tstat = mean(cond.pv);

#--------------------------------------------------------
# pr(p_min \leq M)
pr.pmin.leq.Tstat = length(target.loci)/length(P.min)
# pr(p0 < M) = M
pr.p0.leq.Tstat =  Tstat
# pr(min(p0, p0) \leq M) = 1 - (1 - pr( p_0 \leq M))(1 -pr( p_0 \leq M))
pr.p0min.leq.Tstat = 1 - (1 -   pr.p0.leq.Tstat)*(1 - pr.p0.leq.Tstat )
# pr(p_min \leq M) = 1 - (1 - pr( p_0 \leq M))(1 -pr( p_1 \leq M))
pr.p1.leq.Tstat = max( 1 - (1 -   pr.pmin.leq.Tstat)/(1 - pr.p0.leq.Tstat),
                       pr.p0.leq.Tstat)
# pr(p_1 \leq M) = (1 - pi) pr(p_null \leq M) + pi pr(p_alt \leq M)
pi = (d1 + d2)/p
pr.palt.leq.Tstat = max( (pr.p1.leq.Tstat - (1 - pi) * pr.p0.leq.Tstat)/pi, 1)
#--------------------------------------------------------
# pr(p_0 > p_1 | p_min \leq M) = pr(p_0 > p_1 , p_min \leq M) / pr(p_min <= M)
# pr(p_0 > p_1 , p_min \leq M) \approx pr(p_0 > M , p_1 \leq M)  = pr(p_0 > M ) pr(p_1 \leq M)
pr.p0.geq.p1.given.pmin.leq.Tstat =
  pr.p1.leq.Tstat * (1 - pr.p0.leq.Tstat)/(pr.p1.leq.Tstat + pr.p0.leq.Tstat - pr.p1.leq.Tstat * pr.p0.leq.Tstat )
pr.p0.geq.palt.given.pmin.leq.Tstat =
  pr.palt.leq.Tstat * (1 - pr.p0.leq.Tstat)/(pr.palt.leq.Tstat + pr.p0.leq.Tstat - pr.palt.leq.Tstat * pr.p0.leq.Tstat )
#--------------------------------------------------------
# pr(p_0 <= p_1 | p_min \leq M) = 1 - pr(p_0 > p_1 | p_min \leq M)
pr.p0.leq.p1.given.pmin.leq.Tstat = 1 - pr.p0.geq.p1.given.pmin.leq.Tstat
pr.p0.leq.palt.given.pmin.leq.Tstat = 1 - pr.p0.geq.palt.given.pmin.leq.Tstat
#-----------------------------------------------------------
# pr(p_1 > M, p_min <= M) = pr(p_1 > M, p_0 <= M) = pr(p_1 > M) pr(p_0 <= M)
# pr(p_1 > M | p_min <= M) = pr(p_1 > M, p_min <= M)/pr(p_min <= M)
pr.p1.geq.Tstat.pmin.leq.Tstat = (1 - pr.p1.leq.Tstat)*pr.p0.leq.Tstat
pr.p1.geq.Tstat.given.pmin.leq.Tstat = (1 - pr.p1.leq.Tstat)*pr.p0.leq.Tstat/
  (pr.pmin.leq.Tstat)
pr.palt.geq.Tstat.pmin.leq.Tstat = (1 - pr.palt.leq.Tstat)*pr.p0.leq.Tstat
pr.palt.geq.Tstat.given.pmin.leq.Tstat = (1 - pr.palt.leq.Tstat)*pr.p0.leq.Tstat/
  (pr.pmin.leq.Tstat)

sum.A = length(target.loci) * pr.p0.leq.p1.given.pmin.leq.Tstat -
  sum(1/cond.pv) * pr.p1.geq.Tstat.given.pmin.leq.Tstat

sum.A.2 = (length(target.loci) - length(P.min)* pr.p0min.leq.Tstat)*(
  pr.p0.leq.palt.given.pmin.leq.Tstat -
    pr.p0.geq.Tstat.given.pmin.leq.Tstat * pr.palt.geq.Tstat.given.pmin.leq.Tstat
)

print(c(pr.p0.leq.palt.given.pmin.leq.M,
        pr.p0.geq.M.given.p0.geq.pmin.pmin.leq.M,
        pr.palt.geq.M.given.palt.geq.pmin.pmin.leq.M,
        num_alt))

#return(exp(-max(sum.A.2, 0)))
#return(1-A)
}

dp.to.z= function(mu, px){

  z = qnorm(px/2,lower.tail = F)
  dp = 0.5* dnorm(z, mean = mu)/ dnorm(z, mean = 0) +
    0.5* dnorm(-z, mean = mu)/ dnorm(-z, mean = 0)
  return(dp)

}

pp.to.z= function(mu, px){
  # pr(p < px)
  z = qnorm(px/2,lower.tail = F)

  return(
    1 - (pnorm(z, mean = mu) - pnorm(-z, mean = mu))
  )

}

oracle = function(P1, P2, mu){

  P.min = pmin(P1, P2)
  P.max = pmax(P1, P2)
  M = min(P.max)
  target.loci = which(P.min <= M)
  P.min.target =  P.min[target.loci]
  mu = mu[target.loci]

  pr =sapply(1:length(mu), function(x){

    ((1-M) * dp.to.z(px = P.min.target[x],  mu = mu[x]) +
            1 - pp.to.z(mu = mu[x], px =M )   )/
      ((1-P.min.target[x]) * dp.to.z(px = P.min.target[x],
                                     mu = mu[x]) +
         1 - pp.to.z(mu = mu[x], px =P.min.target[x] ) )


  } )

  pr.approx = (1-M)/(1-P.min.target)
  ratio =pr/ pr.approx
  A = sum(1- ratio)
  prod(ratio)
  exp(-A)
  #stopifnot(!anyNA(pr))
  #pr = ((1-M) * dp.to.z(px = P.min.target,  mu = mu) + 1 - pp.to.z(mu = mu, px =M )   )/
  #  ((1-P.min.target) * dp.to.z(px = P.min.target,  mu = mu) + 1 - pp.to.z(mu = mu, px =P.min.target ) )

  return( (1 - prod(pr)))

}



oracle2 = function(P1, P2, mu){

  P.min = pmin(P1, P2)
  P.max = pmax(P1, P2)
  M = min(P.max)
  target.loci = which(P.min <= M)
  P.min.target =  P.min[target.loci]
  mu = mu[target.loci]

  pr =sapply(1:length(mu), function(x){

    ((1-M) * dp.to.z(px = P.min.target[x],  mu = mu[x]) +
       1 - pp.to.z(mu = mu[x], px =M )   )/
      ((1-P.min.target[x]) * dp.to.z(px = P.min.target[x],
                                     mu = mu[x]) +
         1 - pp.to.z(mu = mu[x], px =P.min.target[x] ) )

  } )

  pr.approx = (1-M)/(1-P.min.target)
  ratio =pr/ pr.approx
  A = sum(1- ratio)
  exp(-A)
  #stopifnot(!anyNA(pr))
  #pr = ((1-M) * dp.to.z(px = P.min.target,  mu = mu) + 1 - pp.to.z(mu = mu, px =M )   )/
  #  ((1-P.min.target) * dp.to.z(px = P.min.target,  mu = mu) + 1 - pp.to.z(mu = mu, px =P.min.target ) )
  return( (1 - prod(pr.approx)* exp(-A) ))

}



oracle3 = function(P1, P2, mu){

  P.min = pmin(P1, P2)
  P.max = pmax(P1, P2)
  M = min(P.max)
  target.loci = which(P.min <= M)
  P.min.target =  P.min[target.loci]
  mu = mu[target.loci]

  pr =sapply(1:length(mu), function(x){

    ((1-M) * dp.to.z(px = P.min.target[x],  mu = mu[x]) +
       1 - pp.to.z(mu = mu[x], px =M )   )/
      ((1-P.min.target[x]) * dp.to.z(px = P.min.target[x],
                                     mu = mu[x]) +
         1 - pp.to.z(mu = mu[x], px =P.min.target[x] ) )

  } )

  pr.p0.leq.p1.given.pmin = (1 - pp.to.z(mu = mu, px =P.min.target ))/
    ((1-P.min.target) * dp.to.z(px = P.min.target,mu = mu) +
       1 - pp.to.z(mu = mu, px =P.min.target ) )

  pr.palt.geq.M.given.palt.geq.pmin = (1 - pp.to.z(mu = mu, px =M ))/
    (1 - pp.to.z(mu = mu, px =P.min.target ) )

  pr.p0.geq.M.given.p0.geq.pmin = (1-M)/(1-P.min.target)

  # underlying truth
  A2 = (1 - pr.palt.geq.M.given.palt.geq.pmin/ pr.p0.geq.M.given.p0.geq.pmin)*
    pr.p0.leq.p1.given.pmin
  sum.A = sum(A2)

  #--------------------------------- estimate with the mean
  pr.p0.leq.palt.given.pmin.leq.M =mean(pr.p0.leq.p1.given.pmin[mu != 0])
  pr.p0.geq.M.given.p0.geq.pmin.pmin.leq.M = mean(pr.p0.geq.M.given.p0.geq.pmin[mu != 0])
  pr.palt.geq.M.given.palt.geq.pmin.pmin.leq.M = mean(pr.palt.geq.M.given.palt.geq.pmin[mu != 0])
  num_alt = sum(mu != 0)

  #---------------------------------
  h = (1 - pr.palt.geq.M.given.palt.geq.pmin.pmin.leq.M/ pr.p0.geq.M.given.p0.geq.pmin.pmin.leq.M ) *
    pr.p0.leq.palt.given.pmin.leq.M * num_alt

  # pr(p_0 > p_a | min(p_0, p_a) \leq M)
  1-(1 - M)*pp.to.z(mu = 1, px =M)/(1 - (1 - M)*(1 -pp.to.z(mu = 1, px =M)) ) # what's this formula? I forgot

  mean(M - P.min.target[mu != 0])

  M *(1 -pp.to.z(mu = 1, px =M))/(M - mean(M - P.min.target) ) # not work
  M *(1 -pp.to.z(mu = 1, px =M))/(M +  (1-M)*pp.to.z(mu = 1, px =M) ) #  work!

  #h = (1 - mean(pr.palt.geq.M.given.palt.geq.pmin[mu != 0])/
  #       mean(pr.p0.geq.M.given.p0.geq.pmin[mu != 0]))*
  #  pr.p0.leq.palt.given.pmin.leq.M * sum(mu != 0)


  #---------
  # The next question is how to estimate these quantities?

  if(sum(mu != 0) == 0){h = 0}

  # need input pi: proportion of non-null
  pi = 0.02

  pr.pmin.leq.Tstat.est = length(target.loci)/length(P1)

  pr.p0.leq.Tstat =  M
  # pr(min(p0, p0) \leq M) = 1 - (1 - pr( p_0 \leq M))(1 -pr( p_0 \leq M))
  pr.p0min.leq.Tstat = 1 - (1 -   pr.p0.leq.Tstat)*(1 - pr.p0.leq.Tstat )
  # pr(p_min \leq M) = 1 - (1 - pr( p_0 \leq M))(1 -pr( p_1 \leq M))
  pr.p1.leq.Tstat.est = max( 1 - (1 -   pr.pmin.leq.Tstat.est)/(1 - pr.p0.leq.Tstat),
                         pr.p0.leq.Tstat)

  pr.p0.geq.p1.given.pmin.leq.Tstat.est =
    pr.p1.leq.Tstat.est * (1 - pr.p0.leq.Tstat)/(pr.p1.leq.Tstat.est + pr.p0.leq.Tstat - pr.p1.leq.Tstat * pr.p0.leq.Tstat )

  num_alt.est =  (length(target.loci) - (1 - pi) * length(P.min)* pr.p0min.leq.Tstat)


  pr.approx = (1-M)/(1-P.min.target)

  #stopifnot(!anyNA(pr))
  #pr = ((1-M) * dp.to.z(px = P.min.target,  mu = mu) + 1 - pp.to.z(mu = mu, px =M )   )/
  #  ((1-P.min.target) * dp.to.z(px = P.min.target,  mu = mu) + 1 - pp.to.z(mu = mu, px =P.min.target ) )
  return( (1 - prod(pr.approx)* exp(-h) ))

}



oracle4 = function(P1, P2){

  P.min = pmin(P1, P2)
  P.max = pmax(P1, P2)
  M = min(P.max)
  target.loci = which(P.min <= M)
  P.min.target =  P.min[target.loci]
  pr.approx = (1-M)/(1-P.min.target)
  target.num = length(target.loci)

  #---------  estimate: need to specify ---------------------
  pi = 0.02
  pr.p0.leq.Tstat =  M
  # pr(min(p0, p0) \leq M) = 1 - (1 - pr( p_0 \leq M))(1 -pr( p_0 \leq M))
  pr.p0min.leq.Tstat = 1 - (1 -   pr.p0.leq.Tstat)*(1 - pr.p0.leq.Tstat )
  num_alt.hat =  (target.num - (1 - pi) * length(P.min)* pr.p0min.leq.Tstat)
  pr.pmin.leq.Tstat.est = target.num/length(P1)
  # pr(p_min \leq M) = 1 - (1 - pr( p_0 \leq M))(1 -pr( p_1 \leq M))
  pr.p1.leq.Tstat.est = max( 1 - (1 -   pr.pmin.leq.Tstat.est)/(1 - pr.p0.leq.Tstat),
                             pr.p0.leq.Tstat)
  pr.palt.leq.Tstat = min( (pr.p1.leq.Tstat.est - (1 - pi) * pr.p0.leq.Tstat)/pi, 1)

  #--------------------------------- estimate
  pr.p0.leq.palt.given.pmin.leq.M.hat = 1-(1 - M)*pr.palt.leq.Tstat/(1 - (1 - M)*(1 -pr.palt.leq.Tstat) )
  pr.p0.geq.M.given.p0.geq.pmin.pmin.leq.M.hat = (sum(pr.approx) - 2*(1-M)/(2-M) *(target.num - num_alt.hat))/num_alt.hat
  pr.palt.geq.M.palt.leq.p0.given.pmin.leq.M.hat = M *(1 -pp.to.z(mu = 1, px =M))/(M +  (1-M)*pr.palt.leq.Tstat )

  #pr.palt.geq.M.given.palt.geq.pmin.pmin.leq.M.hat = mean(pr.palt.geq.M.given.palt.geq.pmin[mu != 0])

  #---------------------------------
  h.hat = (pr.p0.leq.palt.given.pmin.leq.M.hat - pr.palt.geq.M.palt.leq.p0.given.pmin.leq.M.hat/
             pr.p0.geq.M.given.p0.geq.pmin.pmin.leq.M.hat )  * num_alt.hat


  pr.approx = (1-M)/(1-P.min.target)

  #stopifnot(!anyNA(pr))
  #pr = ((1-M) * dp.to.z(px = P.min.target,  mu = mu) + 1 - pp.to.z(mu = mu, px =M )   )/
  #  ((1-P.min.target) * dp.to.z(px = P.min.target,  mu = mu) + 1 - pp.to.z(mu = mu, px =P.min.target ) )
  return( (1 - prod(pr.approx)* exp(-h.hat) ))

}

oracle5 = function(P1, P2){

  P.min = pmin(P1, P2)
  P.max = pmax(P1, P2)
  M = min(P.max)
  target.loci = which(P.min <= M)
  P.min.target =  P.min[target.loci]
  pr.approx = (1-M)/(1-P.min.target)
  target.num = length(target.loci)

  #---------  estimate: need to specify ---------------------
  pi = 0.02
  pr.p0.leq.Tstat =  M
  # pr(min(p0, p0) \leq M) = 1 - (1 - pr( p_0 \leq M))(1 -pr( p_0 \leq M))
  pr.p0min.leq.Tstat = 1 - (1 -   pr.p0.leq.Tstat)*(1 - pr.p0.leq.Tstat )
  num_alt.hat =  (target.num - (1 - pi) * length(P.min)* pr.p0min.leq.Tstat)
  pr.pmin.leq.Tstat.est = target.num/length(P1)
  # pr(p_min \leq M) = 1 - (1 - pr( p_0 \leq M))(1 -pr( p_1 \leq M))
  pr.p1.leq.Tstat.est = max( 1 - (1 -   pr.pmin.leq.Tstat.est)/(1 - pr.p0.leq.Tstat),
                             pr.p0.leq.Tstat)
  pr.palt.leq.Tstat = min( (pr.p1.leq.Tstat.est - (1 - pi) * pr.p0.leq.Tstat)/pi, 1)

  #--------------------------------- estimate
  pr.p0.leq.palt.given.pmin.leq.M.hat = 1-(1 - M)*pr.palt.leq.Tstat/(1 - (1 - M)*(1 -pr.palt.leq.Tstat) )
  pr.p0.geq.M.given.p0.geq.pmin.pmin.leq.M.hat = (sum(pr.approx) - 2*(1-M)/(2-M) *(target.num - num_alt.hat))/num_alt.hat
  pr.palt.geq.M.palt.leq.p0.given.pmin.leq.M.hat = M *(1 -pp.to.z(mu = 1, px =M))/(M +  (1-M)*pr.palt.leq.Tstat )

  #pr.palt.geq.M.given.palt.geq.pmin.pmin.leq.M.hat = mean(pr.palt.geq.M.given.palt.geq.pmin[mu != 0])

  #---------------------------------
  h.hat = (pr.p0.leq.palt.given.pmin.leq.M.hat - pr.palt.geq.M.palt.leq.p0.given.pmin.leq.M.hat/
             pr.p0.geq.M.given.p0.geq.pmin.pmin.leq.M.hat )  * num_alt.hat


  pr.approx = (1-M)/(1-P.min.target)

  #stopifnot(!anyNA(pr))
  #pr = ((1-M) * dp.to.z(px = P.min.target,  mu = mu) + 1 - pp.to.z(mu = mu, px =M )   )/
  #  ((1-P.min.target) * dp.to.z(px = P.min.target,  mu = mu) + 1 - pp.to.z(mu = mu, px =P.min.target ) )
  return( (1 - prod(pr.approx)* exp(-h.hat) ))

}


oracle2.update = function(P1, P2, mu){

  P.min = pmin(P1, P2)
  P.max = pmax(P1, P2)
  M = min(P.max)
  target.loci = which(P.min <= M)
  P.min.target =  P.min[target.loci]
  mu = mu[target.loci]

  I.1 = which(mu != 0)
  I.0 = which(mu == 0)

  pr.approx = (1-M)/(1-P.min.target)
  F_min = pp.to.z(mu = mu, px =M) + M - M*(pp.to.z(mu = mu, px =M))
  A.1 = sum((1 - F_min[I.1]) /F_min[I.1])
  A.2 = -log(1 - M) - (M - M^2/2)/(1 - M)^2
  A.3 = sum( log(1 - P.min.target[I.1]))

  A = A.1*A.2 - A.3

  return( (1 - prod(pr.approx)* exp(-A) ))

}


oracle3.update = function(P1, P2, mu){

  P.min = pmin(P1, P2)
  P.max = pmax(P1, P2)
  M = min(P.max)
  target.loci = which(P.min <= M)
  P.min.target =  P.min[target.loci]
  mu = mu[target.loci]

  I.1 = which(mu != 0)
  I.0 = which(mu == 0)

  pr.approx = (1-M)/(1-P.min.target)
  F_min = pp.to.z(mu = mu, px =M) + M - M*(pp.to.z(mu = mu, px =M))
  A.1 = sum((1 - F_min) /F_min)
  A.2 = -log(1 - M) - (M - M^2/2)/(1 - M)^2
  A.3 = sum(log(1 - P.min.target))

  A = A.1*A.2 - A.3

  return( (1 - prod(pr.approx)* exp(-A) ))

}


oracle3.estimate = function(P1, P2){

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

  return( (1 - prod(pr.approx)* exp(-A) ))

}
cond.prob = function(mu.1, mu.2, cond.p.min){

  cond.z.min = qnorm(cond.p.min/2,lower.tail = F)

  P.min = pmin(P1, P2)

  prob = mean(P1[P.min < thresh] > P2[P.min < thresh])

  return(prob)

}

cond.prob = function(P1, P2, thresh){

  P.min = pmin(P1, P2)

  prob = mean(P1[P.min < thresh] > P2[P.min < thresh])

  return(prob)

}

