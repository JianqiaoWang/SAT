
#' Title
#'
#' @param P1
#' @param P2
#' @param snplist matched SNP name
#' @param CHR matched choromsome name
#' @param POS matched snp position, same length as P1 and P2
#' @param dist_thresh choosed distance threshold
#' @param ref.geno reference genotype to calculate the
#' @param ref.SNP reference genotype column name
#' @param ref.ind.row selected row of reference genotype individuals
#'
#' @return
#' @export
#'
#' @examples
#'
#'
#'
#' MinMaxP.geno.detect(sumstat.comb$P.gwas, sumstat.comb$P.gtex, bdist = "BP", dist_thresh = 500000, ref.geno = geno, ref.geno.ind.col = ref.geno.ind.col)

geno_detect = function(P1, P2, snplist, CHR, POS, dist_thresh = 0.5, ref.geno = NULL,
                               ref.SNP, ref.row = 1:nrow(ref.geno) ){

  # The point is whether the CHR and POS should be aligned with reference genotype data
  # require the matched ID between geno and (P1, P2)
  P.max = pmax(P1, P2); P.min = pmin(P1, P2); M = min(P.max)
  I.m = which(P.max <= M)
  C_0_m = Cmax_null(t = M, P.min = P.min)
  snp_info = Find.block.distance(SNP = snplist[I.m], CHR = CHR[I.m],
                                 POS = POS[I.m], dist_thresh = dist_thresh)

  snp_info$ref_index =   match(snp_info$SNP, ref.SNP) # suppose we have a full match
  Block.ref.list = split(snp_info$ref_index, snp_info$block_id)

  neff_vec = sapply(Block.ref.list, function(x){
    if(length(x) == 1){ return(1)
    }else{
      corr = cor(ref.geno[ref.row,x], use = "pairwise.complete.obs")
      cal_neff(corr, block.thresh = 1 )
    }
  } )
  neff = sum(neff_vec)
  p.vec.dep = 1 - exp(- neff * C_0_m)
  return(list(PV = p.vec.dep, Ne = neff))
}

#' Find the block by the snp distance and position
#' @export
Find.block.distance = function(SNP, CHR, POS = NULL, dist_thresh = 500000){

  #determine number of blocks
  snp_data <- data.frame(
    SNP = SNP,
    CHR = CHR,
    POS = POS
  )
  snp_data= snp_data[order(snp_data$CHR, snp_data$POS),]
  snp_data$POS_diff = snp_data$POS - c(snp_data$POS[1], snp_data$POS[-length(snp_data$POS)] )
  snp_data$CHR_diff = snp_data$CHR - c(snp_data$CHR[1], snp_data$CHR[-length(snp_data$CHR)] )
  snp_data$block_incre = pmax( abs(snp_data$POS_diff) > distance_threshold, abs(snp_data$CHR_diff) != 0 )
  snp_data$block_id = cumsum(snp_data$block_incre)+1

  #determine which snp belonging to the block
  return(snp_data)
}

#' Find the block by the snp distance and position
#' @export
cal_neff = function(corr, block.thresh = 1){
  corr1= r_to_rho(corr)
  corr1[corr1 < 0] = 0
  if(block.thresh < 0.999){
    corr1 = RemoveRedunt(corr1,
                         block.thresh = block.thresh)
  }
  evs = eigen(corr1, only.values = T)$values
  abs.evs <- abs(evs) + sqrt(.Machine$double.eps)
  k <- length(evs)
  n.block.eff <- k - sum( ifelse(abs.evs >= 1, 1, 0) *(abs.evs -1))
  return(n.block.eff)
}
