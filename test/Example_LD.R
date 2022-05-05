p <- 1*10^4; ## total number of pairs
d1 = 100; d2 = 100; d3 = 0;
X <- c(rep(0,p-(d1 + d2 + d3)),rep(1,d1),rep(2,d2),rep(3,d3));
mu.1 = rep(0, p); mu.2 = rep(0, p);
mu.1[X==1|X==3] = 4; mu.2[X==2|X==3] = 5;

#------------------------ add LD matrix into the generation
#source('HD-Dependent-Z.R')
#source('MaxDetect_p.R')
library(dplyr)
library(Matrix)
nblock = 1000 ; Pblock = rep(p/nblock, nblock) #number of SNP
rho.vector = runif(nblock, min = 0.9, max = 0.95)
SIGMA = lapply(1:nblock, function(i) {
  rho.vector[i] ^ abs(outer(1:Pblock[i], 1:Pblock[i], "-"))
})
SIGMA.SQRT = lapply(1:nblock, function(i) {
  eires = eigen(SIGMA[[i]])
  U = eires$vectors
  Lambda = diag(eires$values)
  return(U %*% sqrt(Lambda) %*% U )
}) %>% bdiag_m()
#---------------
SIGMA = SIGMA %>% bdiag_m()
mu.1 = (SIGMA %*% mu.1) %>% as.vector()
mu.2 = (SIGMA %*% mu.2) %>% as.vector()
# -----------------------
p.vec = list()
p.vec.2 = list()
p.vec.3 = list()
for(i in 1:400){
  Z1 <- mu.1 + GenZ(SIGMA.SQRT); Z2 <- mu.2 + GenZ(SIGMA.SQRT);
  T1 = abs(Z1); T2 = abs(Z2);
  P1 = 2* (1 - pnorm( T1)); P2 = 2* (1 - pnorm(T2));
  p.vec[[i]] = MinDetect(P1, P2, K = 400, F.null = NULL)$p;
  b<-Sys.time()
  print(b-a)
}

