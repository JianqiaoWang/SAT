p <- 1*10^5; ## total number of pairs
d1 = 100; d2 = 100; d3 = 0;
X <- c(rep(0,p-(d1 + d2 + d3)),rep(1,d1),rep(2,d2),rep(3,d3));
mu.1 = rep(0, p); mu.2 = rep(0, p);
mu.1[X==1|X==3] = 4; mu.2[X==2|X==3] = 5;


p.vec = list()
p.vec.2 = list()
p.vec.3 = list()
for(i in 1:400){
  Z1 <- mu.1 + rnorm(p,0,1);
  Z2 <- mu.2 + rnorm(p,0,1); 
  T1 = abs(Z1);T2 = abs(Z2);
  P1 = 2* (1 - pnorm( T1)); P2 = 2* (1 - pnorm(T2));
  p.vec[[i]] = MaxDetect.GWAS.Z(P1, P2, K = 500)$p
  #a<-Sys.time()
  p.vec.2[[i]] = ssa::maxtest(T1,T2)$p;
  #b<-Sys.time()
  #print(b-a)
  #print( p.vec[[i]]) ; print(p.vec.2[[i]])
}
print(mean(p.vec <= 0.05));
print(mean(p.vec.2 <= 0.05));

hist(unlist(p.vec));hist(unlist(p.vec.2));