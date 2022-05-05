#' general detection function
#'
#' @param P1
#' @param P2
#'
#' @return
#' @export
#'
#' @examples
#'
#'
#'



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
