#' calculate the effective number of given locations
#'
#' @param Block_list
#' @param ref.geno
#' @param ref.bed
#' @param block.map
#'
#' @return
#' @export
#'
#' @examples
#'
#'
#'
#'
#'
N_eff = function(block.list, ref.geno = NULL, block.thresh = 0.995, ind.row = rows_along(ref.geno) ){

  n.block = length(block.list)

  n.eff = lapply(block.list, function(x){
    if(length(x) == 1){
      n.block.eff = 1
    }else{
      corr = cor(ref.geno[ind.row,x], use = "pairwise.complete.obs")
      corr1= r_to_rho(corr)
      corr1[corr1 < 0] = 0
      corr1 = RemoveRedunt(corr1,
                           block.thresh = block.thresh)
      evs2 = eigen(corr1, only.values = T)$values
      n.block.eff = meff(eigen = evs2, method = "li2012")
    }
    return(data.frame(n.block.eff = n.block.eff) )

  }) %>%  dplyr::bind_rows() %>% colSums()

  return(n.eff)
}

#' @param Sigma.list
#' @param ref.geno
#' @return
#' @export
#'
#' @examples
#'

N_eff_Sigma = function(Sigma.list,
                 block.thresh = 0.995){

  n.block = length(Sigma.list)

  n.eff = lapply(Sigma.list, function(x){
    if(length(x) == 1){
      n.block.eff = 1
    }else{
      corr1= r_to_rho(x)
      corr1[corr1 < 0] = 0
      if(block.thresh < 0.999){
      corr1 = RemoveRedunt(corr1,
                           block.thresh = block.thresh)
      }
      evs2 = eigen(corr1, only.values = T)$values
      n.block.eff = meff(eigen = evs2, method = "li2012")
    }
    return(data.frame(n.block.eff = n.block.eff) )

  }) %>%  dplyr::bind_rows() %>% colSums()

  return(n.eff)
}


#' @export


N_eff_old = function(var.name.target,
                 ref.geno = NULL,
                 ref.bed = NULL,
                 block.map = NULL,
                 var.name = NULL,
                 output.dir = "./temp/",
                 block.thresh = 0.995,
                 K = 1000,
                 exact = F){
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
  return(n.eff)
}


n.eff.vec <- function(ref.geno, block.list, block.thresh){
  n.eff = lapply(block.list, function(x){
    if(length(x) == 1){
      n.eff.liji.2 = 1
    }else{
      corr = cor(ref.geno[,x], use = "pairwise.complete.obs")
      corr1= r_to_rho(corr)
      corr1[corr1 < 0] = 0
      corr1 = RemoveRedunt(corr1)
      #evs = eigen(corr, only.values = T)$values
      evs2 = eigen(corr1, only.values = T)$values
     # n.eff.liji = SAT::meff(eigen = evs, method = "liji")
      n.eff.liji.2 = meff(eigen = evs2, method = "li2012")
      #n.eff.galwey = SAT::meff(eigen = evs, method = "galwey")
    }

    return(data.frame(n.eff.liji.p = n.eff.liji.2) )

  }) %>%  dplyr::bind_rows() %>% colSums()
  return(n.eff)
}

#------------------------- main function ------------------------
r_to_rho = function(x){
  x = x^2
  y = 0.7723 *x^6 - 1.5659 * x^5 + 1.201* x^4 - 0.2355* x^3 + 0.2184*x^2 + 0.6086*x
  return(y)
}

RemoveRedunt = function(corr, block.thresh=0.995){
  arrind = which(abs(corr) >= block.thresh, arr.ind = T)
  arrind = arrind[ arrind[,1] < arrind[,2], ]
  if(length(arrind) > 0){
    remove.ind = RemoveDuplicate(arrind)
    if(length(remove.ind) > 0) {corr = corr[-remove.ind, -remove.ind]}
  }
  return(corr)
}
#
# # neff calculation
# n.eff.calulate <- function(ref.geno, block.list, thresh = 0.995, colsum = T){
#   n.eff = lapply(block.list, function(x){
#     if(length(x) == 1){
#       n.eff.nyholt =  1
#       n.eff.liji = 1
#       n.eff.liji.2 = 1
#       n.eff.gao = 1
#       n.eff.galwey = 1
#     }else{
#
#       corr = cor(ref.geno[,x], use = "pairwise.complete.obs")
#       arrind = which(abs(corr) >= thresh, arr.ind = T)
#       arrind = arrind[ arrind[,1] < arrind[,2], ]
#       if(length(arrind) > 0){
#         remove.ind = RemoveDuplicate(arrind)
#         if(length(remove.ind) > 0) {corr = corr[-remove.ind, -remove.ind]}
#       }
#       corr1= r_to_rho(corr)
#       #corr1 =  poolr::mvnconv(corr, target = "p", cov2cor = TRUE)
#       #evs = eigen(corr, only.values = T)$values
#       evs2 = eigen(corr1, only.values = T)$values
#       #n.eff.nyholt =  SAT::meff(eigen = evs, method = "nyholt")
#       #n.eff.liji = SAT::meff(eigen = evs, method = "liji")
#       n.eff.liji.2 = SAT::meff(eigen = evs2, method = "li2012")
#       #n.eff.gao = SAT::meff(eigen = evs, method = "gao")
#       #n.eff.galwey = SAT::meff(eigen = evs, method = "galwey")
#     }
#
#     n.eff = n.eff.liji.2
#   #   return( data.frame(#n.eff.nyholt = n.eff.nyholt,
#   #                      #n.eff.liji = n.eff.liji,
#   #                      n.eff.liji.p = n.eff.liji.2
#   #                      #n.eff.gao = n.eff.gao,
#   #                      #n.eff.galwey = n.eff.galwey
#   #   ) )
#   #
#   # }) %>%  dplyr::bind_rows()
#   # if(colsum){n.eff = n.eff %>% colSums()}
#   return(n.eff)
# }
#
#
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
#
#
# n.eff.vec.simu <- function(SIGMA, block.list){
#   n.eff = lapply(block.list, function(x){
#     if(length(x) == 1){
#       n.eff.nyholt =  1
#       n.eff.liji = 1
#       n.eff.liji.2 = 1
#       n.eff.gao = 1
#       n.eff.galwey = 1
#     }else{
#       corr = as.matrix(SIGMA[x, x])
#       corr1 =  poolr::mvnconv(corr, target = "p", cov2cor = TRUE)
#       evs = eigen(corr, only.values = T)$values
#       evs2 = eigen(corr1, only.values = T)$values
#       n.eff.nyholt =  SAT::meff(eigen = evs, method = "nyholt")
#       n.eff.liji = SAT::meff(eigen = evs, method = "liji")
#       n.eff.liji.2 = SAT::meff(eigen = evs2, method = "li2012")
#       n.eff.gao = SAT::meff(eigen = evs, method = "gao")
#       n.eff.galwey = SAT::meff(eigen = evs, method = "galwey")
#     }
#     return( data.frame(n.eff.nyholt = n.eff.nyholt,
#                        n.eff.liji = n.eff.liji,
#                        n.eff.liji.p = n.eff.liji.2,
#                        n.eff.gao = n.eff.gao,
#                        n.eff.galwey = n.eff.galwey
#     )
#     )
#
#   }) %>%  dplyr::bind_rows() %>% colSums()
#   return(n.eff)
# }
#
#
#
meff <- function(R, eigen, method, ...) {

  # match 'method' argument
  method <- match.arg(method, c("nyholt", "liji", "gao", "galwey", "li2012"))

  if (missing(eigen)) {

    # check if 'R' is specified
    if (missing(R))
      stop("Argument 'R' must be specified.", call.=FALSE)

    # checks for 'R' argument
    R <- .check.R(R, checksym = TRUE, checkna = TRUE, checkpd = FALSE, nearpd = FALSE, checkcor = TRUE, checkdiag = TRUE, isbase = FALSE)

    # get eigenvalues of 'R' matrix
    evs <- base::eigen(R)$values

  } else {

    # can pass eigenvalues directly to function via 'eigen'

    if (!is.vector(eigen))
      stop("Argument 'eigen' must be a numeric vector.", call.=FALSE)

    evs <- eigen

  }

  # check if there are negative eigenvalues
  if (any(evs < 0))
    warning(paste0("One or more eigenvalues ", ifelse(missing(eigen), "derived from the 'R' matrix ", ""), "are negative."), call.=FALSE)

  if (method == "nyholt") {

    # effective number of tests (based on Nyholt, 2004)
    k <- length(evs)
    m <- 1 + (k - 1) * (1 - var(evs) / k)

  }

  if (method == "li2012"){

    # effective number of tests (based on Li, 2011)
    # adding a small value to the absolute eigenvalues to overcome numerical imprecisions
    abs.evs <- abs(evs) + sqrt(.Machine$double.eps)
    k <- length(evs)
    m <- k - sum( ifelse(abs.evs >= 1, 1, 0) *(abs.evs -1))

  }

  if (method == "liji") {

    # effective number of tests (based on Li & Ji, 2005)
    # adding a small value to the absolute eigenvalues to overcome numerical imprecisions
    abs.evs <- abs(evs) + sqrt(.Machine$double.eps)
    m <- sum(ifelse(abs.evs >= 1, 1, 0) + (abs.evs - floor(abs.evs)))

  }

  if (method == "gao") {

    # effective number of tests (based on Gao, 2008)

    ddd <- list(...)

    # allow user to specify value of C via ... but otherwise use 0.995
    if (!is.null(ddd$C)) {
      C <- ddd$C
    } else {
      C <- 0.995
    }

    if (C < 0 || C >= 1)
      warning("Value of 'C' should be >= 0 and < 1.", call.=FALSE)

    m <- which(cumsum(sort(evs, decreasing = TRUE)) / sum(evs) > C)[1]

  }

  if (method == "galwey") {

    # if there are negative eigenvalues, inform user that they were set to 0
    if (any(evs < 0)) {
      warning(paste0("Negative eigenvalues ", ifelse(missing(eigen), "derived from the 'R' matrix ", ""), "were set to 0."), call.=FALSE)
      evs[evs < 0] <- 0
    }

    # effective number of tests (based on Galwey, 2009)
    m <- sum(sqrt(evs))^2 / sum(evs)

  }

  #m <- ceiling(m)
  # always round down the estimated value
  return(m)

}

#
# n.eff.vec <- function(ref.geno, block.list){
#   n.eff = lapply(block.list, function(x){
#     corr = cor(ref.geno[, x], use = "complete.obs")
#     corr1 =  poolr::mvnconv(corr, target = "p", cov2cor = TRUE)
#     evs = eigen(corr, only.values = T)$values
#     evs.tr =eigen(corr1, only.values = T)$values
#     n.eff.1 =  SAT::meff(eigen = evs, method = "nyholt")
#     n.eff.2 = SAT::meff(eigen = evs, method = "liji")
#     n.eff.3 = SAT::meff(eigen = evs, method = "gao")
#     n.eff.4 = SAT::meff(eigen = evs, method = "galwey")
#     n.eff.tr.1 =  SAT::meff(eigen = evs.tr, method = "nyholt")
#     n.eff.tr.2 = SAT::meff(eigen = evs.tr, method = "liji")
#     n.eff.tr.3 = SAT::meff(eigen = evs.tr, method = "gao")
#     n.eff.tr.4 = SAT::meff(eigen = evs.tr, method = "galwey")
#
#     return( data.frame(n.eff.1 = n.eff.1, n.eff.2 = n.eff.2, n.eff.3 = n.eff.3 ,
#                        n.eff.4 = n.eff.4,
#                        n.eff.tr.1 = n.eff.tr.1, n.eff.tr.2 = n.eff.tr.2, n.eff.tr.3 = n.eff.tr.3,
#                        n.eff.tr.4 = n.eff.tr.4))
#
#   }) %>%  dplyr::bind_rows() %>% colSums()
#
#   return(n.eff)
# }
