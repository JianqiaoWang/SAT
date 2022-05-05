p <- 1*10^5; ## total number of pairs
d1 = 100
d2 = 100
d3 = 0
X <- c(rep(0,p-(d1 + d2 + d3)),rep(1,d1),rep(2,d2),rep(3,d3));
p.vec = list()
p.vec.2 = list()
for(i in 1:200){
  set.seed(i);
  Z1 <- rnorm(p,0,1); Z1[X==1|X==3] <- rnorm(d1+d3,4,1);
  Z2 <- rnorm(p,0,1); Z2[X==2|X==3] <- rnorm(d2+d3,5,1);
  T1 = abs(Z1)
  T2 = abs(Z2)
  p.vec[[i]] = maxtest(T1,T2)$p;
  a<-Sys.time()
  p.vec.2[[i]] = MaxDetect(T1, T2, K = 500)$p
  b<-Sys.time()
  print(b-a)
  print( p.vec[[i]]) ; print(p.vec.2[[i]])
}
print(mean(p.vec < 0.05));
print(mean(p.vec.2 < 0.05));
hist(unlist(p.vec));hist(unlist(p.vec.2));