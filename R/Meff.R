#' Detect function for the genotype data
#'
#' @param R
#' @export
#'
#' @examples
#'
#'
#'
#'
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
