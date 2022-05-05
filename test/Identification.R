library(dplyr)
library(parallel)
library(doParallel)
devtools::load_all('../../SAT/')
p <- 3*10^5; ## total number of pairs
#d <-as.numeric(commandArgs(trailingOnly = TRUE)[2])
#r <-as.numeric(commandArgs(trailingOnly = TRUE)[3])
#H <-(commandArgs(trailingOnly = TRUE)[4])
d = 1000;r = 0
d1 = d; d2 = d;
a = r;  b = r
set.seed(123)
d3 = 30;
X <- c(rep(0,p-(d1 + d2 + d3)),rep(1,d1),rep(2,d2),rep(3,d3));
mu.1 = rep(0, p); mu.2 = rep(0, p);
mu.1[X==1] = rnorm(d1, a, 1);
mu.2[X==2] =rnorm(d2, b, 1);
mu.1[X==3] = rnorm(d3, 5, 1);
mu.2[X==3] = 0.5 * mu.1[X==3]

source('CompareFunc.R')

sim = function(i){
  set.seed(i)
  Z1 <- mu.1 + rnorm(p,0,1);
  Z2 <- mu.2 + rnorm(p,0,1);
  T1 = abs(Z1);T2 = abs(Z2);
  P1 = 2* (1 - pnorm( T1));
  P2 = 2* (1 - pnorm(T2));
  p.vec = SAT::MinID(P1, P2)
  print(c( p.vec$pv.sparse, p.vec$pv.dense ))
  return(c(p.vec$pv.sparse, p.vec$pv.dense, cmpr.vec))
}
