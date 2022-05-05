#' Detect function
#'
#' @param P1
#' @param P2
#' @param K
#' @param ref.geno
#'
#' @return
#' @export
#'
#' @examples
MinDetect_LD_cor = function(P1, P2, K = 1000, ref.geno = NULL, block.map = NULL, LD.thresh, exact = F){
  # the input is p-value sequence
  # When LD.thresh is 1, it degrades to the independence case
  # When LD.thresh is 0, it degrades to the total dependence case
  vec = vector()
  P.min = pmin(P1, P2)
  P.max = pmax(P1, P2)
  Tstat = min(P.max)
  target.loci = which(P.min <= Tstat)
  P.min.target =  P.min[target.loci]
  ref.geno = ref.geno[,target.loci]

 #--------------------- define the block list -------------------
 if(is.null(block.map) ){
  adj.cor.mat = sapply(1:(ncol(ref.geno)-1), function(x){
    cor(ref.geno[,x], ref.geno[,x+1], use = "pairwise.complete.obs")
  } )
   block.list =  GenerateBlock(adj.cor.mat, rho = LD.thresh)
 }else{

   block.list =  split(1:length(target.loci), block.map[target.loci])

   }
#-----------------------------------------------------------------
   #block.list.05 =  GenerateBlock(adj.cor.mat, rho = LD.thresh)
   #block.list.07 =  GenerateBlock(adj.cor.mat, rho = LD.thresh)

   pr = sapply(block.list, function(x){ # generate block list  probability
     if(length(x) == 1){
       pr = (1 - Tstat)/(1 - P.min.target[x]); # pr( P.max > Tstat |P.min)
     }else{
       pr = pv.block(Tstat,
                     P.min = P.min.target[x],
                     ref.geno = ref.geno[,x], exact = exact) # pr( P.max > Tstat |P.min) within a block
     }
     return(pr)
   })  # pr( P.max > Tstat |P.min)

   p.value = 1 - prod(pr); # calculate the p value for the observed stat

  return(list(M = Tstat, p = p.value))
}

GenerateBlock = function(adj.cor.mat, rho){
  cut.point = which(abs(adj.cor.mat) <= rho) # Identify non-LD block
  if(length(cut.point) == 0){
    block.division = c(0, length(adj.cor.mat)+1)
  }else{
    block.division = c(0, cut.point, length(adj.cor.mat)+1)
  }
  block.list = lapply(1:(length(block.division)-1), function(x){
    return(c( (block.division[x]+1):block.division[x+1]) )
  })
}



