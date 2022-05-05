#' Detect function for the genotype data
#'
#' @param P1
#' @param P2
#' @param K
#' @param F.null
#'
#' @return
#' @export
#'
#' @examples
#'
#'
#'

MinMaxP.geno = function(P1, P2,
                   ref.geno = NULL,
                   ref.bed = NULL,
                   block.map = NULL,
                   var.name = NULL,
                   output.dir = "./temp/",
                   block.thresh = 0.995,
                   K = 1000,
                   exact = F){
  # the input is p-value sequence
  # geno is a bigsnpr object without missing values

  P.min = pmin(P1, P2)
  P.max = pmax(P1, P2)
  M = min(P.max)
  Tstat = M
  F_min_M = mean(P.min <= M)
  I.m = which(P.min <= M)
  I.m.size = length(I.m)
  P.min.target =  P.min[I.m]
  if(is.null(var.name)){
    var.name.target = I.m
  }else{
    var.name.target =  var.name[I.m]
  }
  names(P.min.target) = var.name.target
  pr.sparse = pr.sparse.ind(M, P.min.target)
  #pr.dense = rep(pr.dense.ind(M, F_min_M), I.m.size)
  log.pr.dense = log(1 - M)/F_min_M + ((1- F_min_M)/F_min_M) * (M - M^2/2)/(1 - M)^2

  if(is.null(ref.geno) & is.null(ref.bed) ){
    cat("Ignoring the dependence")
    return(list(M = M, pv.sparse = 1 - prod(pr.sparse),
                pv.dense = 1 - exp(I.m.size * log.pr.dense) ))
  }
  #------------------ read the genotype file ------------------------
  if(!is.null(ref.geno)){
    loc.target = match(var.name.target, ref.geno$map$marker.ID)
    ref.geno = ref.geno$genotypes[,loc.target] %>% as.matrix
    colnames(ref.geno) = var.name.target
  }else{

  if(!is.null(ref.bed)){
     ref.geno = snpStats::read.plink( bed = ref.bed,
                                      select.snps = var.name.target )
     ref.geno = as(ref.geno$genotypes, "numeric")
    }
  }
  #--------------------- define the block list -------------------
  if(is.null(block.map)){
  block.list = Find.plink.block(snplist = var.name.target,
                                ref.bed = ref.bed,
                                output.dir = output.dir)
  }else{
    block.list = split(var.name.target, block.map[I.m])
  }
  n.block = length(block.list)
  n.eff = n.eff.vec(ref.geno, block.list, block.thresh = block.thresh)
  p.value = 1 - exp(n.eff * log.pr.dense)
  return(list(M = Tstat, p = p.value, pr.dense = exp(I.m.size * log.pr.dense), I.m.size = (I.m.size),
              n.block = n.block, n.eff = n.eff))
}

Adj.func = function(P1, P2){
  P.min = pmin(P1, P2)
  P.max = pmax(P1, P2)
  M = min(P.max)
  I.m = which(P.min <= M)
  P.min.target =  P.min[I.m]
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

  return( exp(-A) )

}

n.eff.vec <- function(ref.geno, block.list, block.thresh){
  n.eff = lapply(block.list, function(x){
    if(length(x) == 1){
      n.eff.liji = 1
      n.eff.liji.2 = 1
      n.eff.galwey = 1
    }else{
      corr = cor(ref.geno[,x], use = "pairwise.complete.obs")
      corr1= r_to_rho(corr)
      corr1[corr1 < 0] = 0
      corr1 = RemoveRedunt(corr1)
      evs = eigen(corr, only.values = T)$values
      evs2 = eigen(corr1, only.values = T)$values
      n.eff.liji = SAT::meff(eigen = evs, method = "liji")
      n.eff.liji.2 = SAT::meff(eigen = evs2, method = "li2012")
      n.eff.galwey = SAT::meff(eigen = evs, method = "galwey")
    }

    return( data.frame(n.eff.liji = n.eff.liji,
                       n.eff.liji.p = n.eff.liji.2,
                       n.eff.galwey = n.eff.galwey
    ) )

  }) %>%  dplyr::bind_rows() %>% colSums()
  return(n.eff)
}

#------------------------- main function ------------------------

r_to_rho = function(x){
  x = x^2
  y = 0.7723 *x^6 - 1.5659 * x^5 + 1.201* x^4 - 0.2355* x^3 + 0.2184*x^2 + 0.6086*x
  return(y)
}

RemoveRedunt = function(corr, block.thresh=0.98){
  arrind = which(abs(corr) >= block.thresh, arr.ind = T)
  arrind = arrind[ arrind[,1] < arrind[,2], ]
  if(length(arrind) > 0){
    remove.ind = RemoveDuplicate(arrind)
    if(length(remove.ind) > 0) {corr = corr[-remove.ind, -remove.ind]}
  }
  return(corr)
}

n.eff.calulate <- function(ref.geno, block.list, thresh = 0.995, colsum = T){
  n.eff = lapply(block.list, function(x){
    if(length(x) == 1){
      n.eff.nyholt =  1
      n.eff.liji = 1
      n.eff.liji.2 = 1
      n.eff.gao = 1
      n.eff.galwey = 1
    }else{
      corr = cor(ref.geno[,x], use = "pairwise.complete.obs")
      arrind = which(abs(corr) >= thresh, arr.ind = T)
      arrind = arrind[ arrind[,1] < arrind[,2], ]
      if(length(arrind) > 0){
        remove.ind = RemoveDuplicate(arrind)
        if(length(remove.ind) > 0) {corr = corr[-remove.ind, -remove.ind]}
      }
      corr1= r_to_rho(corr)
      #corr1 =  poolr::mvnconv(corr, target = "p", cov2cor = TRUE)
      evs = eigen(corr, only.values = T)$values
      evs2 = eigen(corr1, only.values = T)$values
      n.eff.nyholt =  SAT::meff(eigen = evs, method = "nyholt")
      n.eff.liji = SAT::meff(eigen = evs, method = "liji")
      n.eff.liji.2 = SAT::meff(eigen = evs2, method = "li2012")
      n.eff.gao = SAT::meff(eigen = evs, method = "gao")
      n.eff.galwey = SAT::meff(eigen = evs, method = "galwey")
    }

    return( data.frame(n.eff.nyholt = n.eff.nyholt,
                       n.eff.liji = n.eff.liji,
                       n.eff.liji.p = n.eff.liji.2,
                       n.eff.gao = n.eff.gao,
                       n.eff.galwey = n.eff.galwey
    ) )

  }) %>%  dplyr::bind_rows()
  if(colsum){n.eff = n.eff %>% colSums()}
  return(n.eff)
}


RemoveDuplicate = function(index.arr){

  Final.Remove.Ind = vector()
  # index.arr = index.arr[index.arr[,1] < index.arr[,2],] # consider the (i, j) with i < j
  # if(length(index.arr) ==0){ return(Final.Remove.Ind) }
  if(length(index.arr)==2) {
    Final.Remove.Ind = index.arr[2]
    return(Final.Remove.Ind)
  }
  index.arr = index.arr[order(index.arr[,1],decreasing=FALSE),] # first sort the index arr by first colums
  while(1)
  {
    if(length(index.arr)==2) {
      Final.Remove.Ind = c(Final.Remove.Ind, index.arr[2])
      break
    }
    if(length(index.arr) ==0){break}
    b = index.arr[1,1]
    Final.Remove.Ind = c(Final.Remove.Ind, index.arr[(index.arr[,1] == b),2]) # take the row which first element == b
    remove.index.row.logical = index.arr[,1] %in% Final.Remove.Ind | index.arr[,2] %in% Final.Remove.Ind
    index.arr = index.arr[!remove.index.row.logical, ] #
  }
  return(Final.Remove.Ind)
}


MinMaxP.simu = function(P1, P2, var.name = NULL,
                        block.map = NULL,
                        corr = NULL){

  # the input is p-value sequence
  P.min = pmin(P1, P2)
  P.max = pmax(P1, P2)
  M = min(P.max)
  Tstat = M
  F_min_M = mean(P.min <= M)
  I.m = which(P.min <= M)
  I.m.size = length(I.m)
  P.min.target =  P.min[I.m]
  if(is.null(var.name)){
    var.name.target = I.m
  }else{
    var.name.target =  var.name[I.m]
  }
  names(P.min.target) = var.name.target
  pr.sparse = pr.sparse.ind(M, P.min.target)
  pr.dense = rep(pr.dense.ind(M, F_min_M), I.m.size)

  #Adj = Adj.func(P1, P2)
  #A = -log(Adj)

  if(is.null(block.map) & is.null(corr) ){
    cat("Ignoring the dependence")
    return(list(M = M, pv.sparse = 1 - prod(pr.sparse),
                pv.dense = 1 - prod(pr.dense)))
  }
  #------------------ read the genotype file ------------------------
  block.list = split(var.name.target, block.map[I.m])
  n.block = length(block.list)
  #--------------------------------------------------------------
  n.eff = n.eff.vec.simu(SIGMA = corr, block.list = block.list)
  p.value = 1 - ( prod(pr.dense) )^(n.eff/I.m.size)
  return(list(M = Tstat, p = p.value, pr.dense = prod(pr.dense),  I.m.len = length(I.m),
              n.block = n.block, n.eff = n.eff))
}


n.eff.vec.simu <- function(SIGMA, block.list){
  n.eff = lapply(block.list, function(x){
    if(length(x) == 1){
      n.eff.nyholt =  1
      n.eff.liji = 1
      n.eff.liji.2 = 1
      n.eff.gao = 1
      n.eff.galwey = 1
    }else{
    corr = as.matrix(SIGMA[x, x])
    corr1 =  poolr::mvnconv(corr, target = "p", cov2cor = TRUE)
    evs = eigen(corr, only.values = T)$values
    evs2 = eigen(corr1, only.values = T)$values
    n.eff.nyholt =  SAT::meff(eigen = evs, method = "nyholt")
    n.eff.liji = SAT::meff(eigen = evs, method = "liji")
    n.eff.liji.2 = SAT::meff(eigen = evs2, method = "li2012")
    n.eff.gao = SAT::meff(eigen = evs, method = "gao")
    n.eff.galwey = SAT::meff(eigen = evs, method = "galwey")
    }
    return( data.frame(n.eff.nyholt = n.eff.nyholt,
                       n.eff.liji = n.eff.liji,
                       n.eff.liji.p = n.eff.liji.2,
                       n.eff.gao = n.eff.gao,
                       n.eff.galwey = n.eff.galwey
                       )
            )

  }) %>%  dplyr::bind_rows() %>% colSums()
  return(n.eff)
}


