CF = function(x){ 2*pnorm(x) -1}
FDR = function(T1, T2, alpha, CF){
  vec = vector()
  T.min = pmin(T1, T2)
  T.max = pmax(T1, T2)
  n= length(T.min)
  t.choice = sort(unique(c( T.max, T.min) ), decreasing = TRUE)[1:1000]
  G = function(t) {sum( T.min > t )/n }
  G <- Vectorize(G, SIMPLIFY = T)
  S = function(t) sum( (T.max > t) * (1 - CF(t))  )/n
  #S = function(t) mean( (T1 > t)) * mean( (T2 > t))
  S <- Vectorize(S, SIMPLIFY = T)
  FDR.t =  S(t.choice)/ max(G(t.choice) , 1/length(T.max))
  # Take t as the first such that FDR.t < alpha
  if(sum( FDR.t <= alpha ) > 0){
    t_opt = min(t.choice[FDR.t <= alpha])
    rejected.test = which(T.min >= t_opt)
  }else{
    return(0)
  }
}