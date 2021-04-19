p <- 1*10^5; ## total number of pairs
d1 = 100; d2 = 100; d3 = 0;
X <- c(rep(0,p-(d1 + d2 + d3)),rep(1,d1),rep(2,d2),rep(3,d3));
mu.1 = rep(0, p); mu.2 = rep(0, p);
mu.1[X==1|X==3] = 4; mu.2[X==2|X==3] = 5;
res = vector()
for(d1 in c(10, 50, 100) ){
  for(d2 in c(10, 50, 100)){
    for( a  in c(3,4) ){
      for(b in c( 4,5 )){

mu.1[X==1|X==3] = rnorm(d1+d3, a, 1)
mu.2[X==2|X==3] = rnorm(d2+d3, b, 1)
p.vec = list()
p.vec.2 = list()
p.vec.3 = list()
for(i in 1:600){
  Z1 <- mu.1 + rnorm(p,0,1);
  Z2 <- mu.2 + rnorm(p,0,1);
  T1 = abs(Z1);T2 = abs(Z2);
  P1 = 2* (1 - pnorm( T1)); P2 = 2* (1 - pnorm(T2));
  p.vec[[i]] = SAT::MinDetect(P1, P2)$p
  p.vec.2[[i]] = ssa::maxtest(T1,T2)$p;
}
print(mean(p.vec <= 0.05));
print(mean(p.vec.2 <= 0.05));
res = rbind(res, c(d1, d2, a, b, mean(p.vec <= 0.05), mean(p.vec.2 <= 0.05) ) )

}}}}
res = as.data.frame(res)
colnames(res) = c("d1", "d2", "a", "b", "SAT", "SSA")
save(res, file = "res.simu1.null.Rdata")
