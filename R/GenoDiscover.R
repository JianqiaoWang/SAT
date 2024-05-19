
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

geno_discover = function(P1, P2, snplist, CHR, POS, dist_thresh = 0.5, ref.geno = NULL,
                         ref.SNP, ref.row = 1:nrow(ref.geno) ){

  # Initial cut-off value calculation and candidate set
  tau = SAT::MinMaxP.discov(P.max = P.max, P.min = P.min, method = method, alpha = alpha)$tau
  if(is.null(tau)){ tau = min(P.max)}
  I.m = which(P.min <= tau)
  snp_info = Find.block.distance(SNP = snplist[I.m], CHR = CHR[I.m], POS = POS[I.m], dist_thresh = dist_thresh)
  snp_info$Im = I.m

  # it has to match the index of ref_geno
  snp_info$ref_index =   match(snp_info$SNP, ref.SNP) # suppose we have a full match
  snp_info$eff_weight = rep(1, nrow(snp_info))


  # assign the initial effective weight, where the sum is effective number
  Block.ref.list = split(snp_info$ref_index, snp_info$block_id)
  neff_vec_list = N_eff_list(block.list = Block.ref.list, ref.geno  = ref.geno,
                             block.thresh = 1, ind.row =  ref.row)
  eff_weight = unlist(neff_vec_list)/lengths(Block.ref.list)
  if(length( Block.ref.list ) != 0){
    snp_info$eff_weight = rep(eff_weight, times = lengths(Block.ref.list))
  }

  # if tau only gives single discover set, check whether the  p-value passes
  if(sum(P.max <= tau) == 1){
    C_0_m = Cmax_null(t = tau, P.min = P.min)
    neff = sum(snp_info$eff_weight)
    Pv = 1 - exp(- neff * C_0_m)
    if(Pv > alpha){ return(list(S = NULL, tau = NULL))}
  }

  # iteration
  eff_weight_all = rep(1, length(P.max))
  eff_weight_all[I.m] = snp_info$eff_weight
  tau = MinMaxP.discov.FDR.weight(P.max = P.max, P.min = P.min, weight = eff_weight_all,  alpha = alpha )$tau
  S = which(P.max <= tau)

  return(list(S = S, tau = tau))
}
