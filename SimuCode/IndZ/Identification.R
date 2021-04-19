p <- 1*10^5; ## total number of pairs
# d1 = 100; d2 = 100; d3 = 1;
# X <- c(rep(0,p-(d1 + d2 + d3)),rep(1,d1),rep(2,d2),rep(3,d3));
mu.1 = rep(0, p); mu.2 = rep(0, p);
#mu.1[X==1|X==3] = 4; mu.2[X==2|X==3] = 5;
res = vector()
rep = 20

bunch_of_tests <- function(S, I.null){
  #family-wise error
  FWE <- any(S %in% I.null)
  # false discovery proportion
  FDP <- sum(S %in% I.null)/max(1, length(S))
  # true discovery proportion
  TPP <- sum(!(S %in% I.null))/max(1, p -length(I.null))
  return(c(FWE, FDP, TPP))
}

for(d1 in c(50, 100) ){
  for(d2 in c(50, 100)){
    for(d3 in c(10, 50, 100)){
    for( a  in c(4) ){
      for(b in c(4 )){
        mu.1 = rep(0, p); mu.2 = rep(0, p);
      X <- c(rep(0,p-(d1 + d2 + d3)),
             rep(1,d1),rep(2,d2),rep(3,d3));
      I.null = which(X != 3)
      mu.1[X==1|X==3] = rnorm(d1+d3, a, 1);
      mu.2[X==2|X==3] = rnorm(d2+d3, b, 1);
      sum.vec = list(); sum.vec.2 = list();
      sum.vec.3 = list()
      for(i in 1:rep){
        Z1 <- mu.1 + rnorm(p,0,1); Z2 <- mu.2 + rnorm(p,0,1);
        T1 = abs(Z1);T2 = abs(Z2);
        P1 = 2* (1 - pnorm( T1)); P2 = 2* (1 - pnorm(T2));
        S.Reject = SAT::MinID(P1, P2, type = "FWER", alpha = 0.05)$S
        sum.vec[[i]] = bunch_of_tests(S = S.Reject, I.null = I.null)
        temp = vector()
        for(thresh in c( 10^(-2), 10^(-3), 10^(-4), 10^(-5) )){
        S.Reject = which(P1 < thresh & P2 < thresh)
        temp = c(temp, bunch_of_tests(S = S.Reject, I.null = I.null) )
        }
        sum.vec.2[[i]] = temp
     }
      res.1 = colMeans(do.call(rbind, sum.vec) )
      res.2 = colMeans(do.call(rbind, sum.vec.2) )
      res = rbind(res, c(d1, d2,  d3, a, b,
                         res.1, res.2 ) )

}}}}}
res = as.data.frame(res)
colnames(res) = c("d1", "d2", "d3", "a", "b",
                  paste0(c("FWE", "FDP", "TPP"), ".SAT" ),
                  paste0(c("FWE", "FDP", "TPP"), -2),
                  paste0(c("FWE", "FDP", "TPP"), -3),
                  paste0(c("FWE", "FDP", "TPP"), -4),
                  paste0(c("FWE", "FDP", "TPP"), -5))
save(res, file = "simu1.id.FWER.Rdata")
res[,-c(1:5)] = 100 * res[,-c(1:5)]
knitr::kable(res, format = "latex", digits = 1, booktab = T)
