p <- 100000; ## total number of pairs
s = 20
X <- c(rep(0,p-3*s),rep(1,s),rep(2,s),rep(3,s));
## X=0: no signal in either sequence of tests
## X=1: signal in sequence 1 only
## X=2: signal in sequence 2 only
## X=3: simultaneous signal

TestStat = function(s1, t1){
  
  T.stat = max(pmin(s1,t1))
  
}

Permute = function(T1,T2,thresh, Tstat,  K = 300){
  vec = vector()
  for(i in 1:K){
  T1.permute = T1[sample(length(T1),replace = F)]
  
  adj.fac = thresh/pmax(T1.permute + T2, thresh) 
  
  vec = c(vec,
          max(pmin(T1.permute ,T2)*adj.fac)
  )
  }
  p.value = sum(vec >= Tstat)/length(vec)
  return(p.value)
  }

Permute3 = function(T1, T2, Tstat,K = 500){
  vec = vector()
  for(i in 1:K){
   # T1.0.min = pmin(T1[sample(length(T1),replace = F)],
    #                abs(rnorm(p,0,1)))
   # T2.0.min = pmin(T2[sample(length(T1),replace = F)], 
    #                abs(rnorm(p,0,1)))
    
     T1.0.min = pmin(T1,
                    abs(rnorm(p,0,1)))
     T2.0.min = pmin(T2, 
                    abs(rnorm(p,0,1)))
    pair.min = c(T1.0.min,
                 T2.0.min)
   # pair.sum = T1.permute^2 + T2^2
    ind = which.max(pair.min)
    vec = c(vec, pair.min[ind] )
  }
  p.value = sum(vec >= Tstat)/length(vec)
  return(p.value)
}



Permute4 = function(T1, T2, Tstat,K = 200){
  vec = vector()
  for(i in 1:K){
    Z = abs(rnorm(p,0,1))
    T1.0.min = pmin(T1,
                    Z)
    T2.0.min = pmin(T2, Z)
    pair.min = pmax(T1.0.min,
                 T2.0.min)
    ind = which.max(pair.min)
    vec = c(vec, pair.min[ind] )
  }
  p.value = sum(vec >= Tstat)/length(vec)
  return(p.value)
}




Permute5 = function(T1, T2, Tstat,K = 200){
  vec = vector()
  for(i in 1:K){
    # T1.0.min = pmin(T1[sample(length(T1),replace = F)],
    #                abs(rnorm(p,0,1)))
    # T2.0.min = pmin(T2[sample(length(T1),replace = F)], 
    #                abs(rnorm(p,0,1)))
    
    Z = pmin(T1, T2)[sample(length(T1))]
    T1.0.min = pmin(T1,Z)
    T2.0.min = pmin(T2, Z)
    pair.min = pmax(T1.0.min,
                    T2.0.min)
    # pair.sum = T1.permute^2 + T2^2
    ind = which.max(pair.min)
    vec = c(vec, pair.min[ind])
  }
  p.value = sum(vec >= Tstat)/length(vec)
  return(p.value)
}

Permute2 = function(T1, T2, oracle.X,  K = 300){
  vec = vector()
  while(1){
    g = sample(length(T1))
    if(oracle.X[g] == X[g]){
      next
    }
    T1.permute = T1[sample(length(T1),replace = F)]
    pair.min = pmin(T1.permute^2 ,T2^2)
    pair.sum = T1.permute^2 + T2^2
    ind = which.max(pair.min)
    vec = rbind(vec, c(pair.min[ind], pair.sum[ind]))
  }
  return(vec)
}

method1 = function(){
p.vec = list()
p.vec.2 = list()
for(i in 1:200){
#set.seed(i);
Z1 <- rnorm(p,0,1); Z1[X==1|X==3] <- rnorm(s,6,1);
Z2 <- rnorm(p,0,1); Z2[X==2|X==3] <- rnorm(s,5,1);
p.vec[[i]] = maxtest(abs(Z1),abs(Z2))$p;
#p.vec.2[[i]] = Permute3(abs(Z1),abs(Z2),Tstat = maxtest(abs(Z1),abs(Z2))$M)
p.vec.2[[i]] = Permute4(abs(Z1),abs(Z2),Tstat = maxtest(abs(Z1),abs(Z2))$M)
print(p.vec[[i]]);print(p.vec.2[[i]]);
}
mean(p.vec < 0.05);mean(p.vec.2 < 0.05);
print(mean(p.vec < 0.05));
print(mean(p.vec.2 < 0.05));
hist(unlist(p.vec));hist(unlist(p.vec.2));
}

#method1()

method2 = function(){
  p.vec = list()
  p.vec.2 = list()
  for(i in 1:200){
    #set.seed(i);
    Z1 <- rnorm(p,0,1); Z1[X==1|X==3] <- rnorm(s,3,1);
    Z2 <- rnorm(p,0,1); Z2[X==2|X==3] <- rnorm(s,4,1);
    Z3<-rnorm(p,0,1);
    p.vec[[i]] = maxtest(abs(Z1),abs(Z2))$p;
    #p.vec.2[[i]] = Permute3(abs(Z1),abs(Z2),Tstat = maxtest(abs(Z1),abs(Z2))$M)
    p.vec.2[[i]] = Permute4(abs(Z1),abs(Z2),Tstat = maxtest(abs(Z1),abs(Z2))$M)
    print(p.vec[[i]]);print(p.vec.2[[i]]);
  }
  print(mean(p.vec < 0.05));
  print(mean(p.vec.2 < 0.05));
  hist(unlist(p.vec));hist(unlist(p.vec.2));
}


#method3 = function(){
  p.vec = list()
  p.vec.2 = list()
  for(i in 1:200){
    set.seed(i);
    Z1 <- rnorm(p,0,1); Z1[X==1|X==3] <- rnorm(2*s,4,1);
    Z2 <- rnorm(p,0,1); Z2[X==2|X==3] <- rnorm(2*s,5,1);
    T1 = abs(Z1)
    T2 = abs(Z2)
    reject2 = ssa::nfsdr(cbind(T1, T2), 
                         alpha = 0.05, rescale = F)
    reject =  FDR(T1, T2, 0.05, CF)
    result = c( 1 - mean(reject %in% (p - s):p ),
    sum(reject %in% (p - s):p),
    1 - mean(reject2 %in% (p - s):p ),
    sum(reject2 %in% (p - s):p)
    )
    print(result)
    #reject2 = ssa::nfsdr2_all(T1, T2, alpha = 0.05)
  
  }
  print(mean(p.vec < 0.05));
  print(mean(p.vec.2 < 0.05));
  hist(unlist(p.vec));hist(unlist(p.vec.2));
#}
#vec1 = Permute2(T1, T2, K = 800)
#vec2 = Permute2(abs(rnorm(p,0,1)), T2, K = 800)
#vec3 = Permute2(T1,abs(rnorm(p,0,1)), K = 800)
#hist(vec1[,2], breaks = 50)
#hist(vec2[,2], breaks = 30)
#hist(vec3[,2], breaks = 30)
