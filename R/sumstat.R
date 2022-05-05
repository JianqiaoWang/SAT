
#

clean_header = function(names){
  return(
    names %>% toupper %>% gsub("\\.", "_", .) %>% gsub("\n", "", .)
  )
}

#
# h = reticulate::py_to_r(
# reticulate::py_run_string("
# default_cnames = {
#
#     # RS NUMBER
#     'SNP': 'SNP',
#     'MARKERNAME': 'SNP',
#     'SNPID': 'SNP',
#     'RS': 'SNP',
#     'RSID': 'SNP',
#     'RS_NUMBER': 'SNP',
#     'RS_NUMBERS': 'SNP',
#     # NUMBER OF STUDIES
#     'NSTUDY': 'NSTUDY',
#     'N_STUDY': 'NSTUDY',
#     'NSTUDIES': 'NSTUDY',
#     'N_STUDIES': 'NSTUDY',
#     # P-VALUE
#     'P': 'P',
#     'PVALUE': 'P',
#     'P_VALUE':  'P',
#     'PVAL': 'P',
#     'P_VAL': 'P',
#     'P-VALUE': 'P',
#     'GC_PVALUE': 'P',
#     # ALLELE 1
#     'A1': 'A1',
#     'ALLELE1': 'A1',
#     'ALLELE_1': 'A1',
#     'EFFECT_ALLELE': 'A1',
#     'REFERENCE_ALLELE': 'A1',
#     'INC_ALLELE': 'A1',
#     'EA': 'A1',
#     # ALLELE 2
#     'A2': 'A2',
#     'ALLELE2': 'A2',
#     'ALLELE_2': 'A2',
#     'OTHER_ALLELE': 'A2',
#     'NON_EFFECT_ALLELE': 'A2',
#     'DEC_ALLELE': 'A2',
#     'NEA': 'A2',
#     # N
#     'N': 'N',
#     'NCASE': 'N_CAS',
#     'CASES_N': 'N_CAS',
#     'N_CASE': 'N_CAS',
#     'N_CASES': 'N_CAS',
#     'N_CONTROLS': 'N_CON',
#     'N_CAS': 'N_CAS',
#     'N_CON': 'N_CON',
#     'N_CASE': 'N_CAS',
#     'NCONTROL': 'N_CON',
#     'CONTROLS_N': 'N_CON',
#     'N_CONTROL': 'N_CON',
#     'WEIGHT': 'N',  # metal does this. possibly risky.
#     # SIGNED STATISTICS
#     'ZSCORE': 'Z',
#     'Z-SCORE': 'Z',
#     'GC_ZSCORE': 'Z',
#     'Z': 'Z',
#     'OR': 'OR',
#     'B': 'BETA',
#     'BETA': 'BETA',
#     'STDERR': 'SE',
#     'LOG_ODDS': 'LOG_ODDS',
#     'EFFECTS': 'BETA',
#     'EFFECT': 'BETA',
#     'SIGNED_SUMSTAT': 'SIGNED_SUMSTAT',
#     # INFO
#     'INFO': 'INFO',
#     # MAF
#     'EAF': 'FRQ',
#     'FRQ': 'FRQ',
#     'MAF': 'FRQ',
#     'FRQ_U': 'FRQ',
#     'F_U': 'FRQ',
# }
#
# describe_cname = {
#     'SNP': 'Variant ID (e.g., rs number)',
#     'P': 'p-Value',
#     'A1': 'Allele 1, interpreted as ref allele for signed sumstat.',
#     'A2': 'Allele 2, interpreted as non-ref allele for signed sumstat.',
#     'N': 'Sample size',
#     'N_CAS': 'Number of cases',
#     'N_CON': 'Number of controls',
#     'Z': 'Z-score (0 --> no effect; above 0 --> A1 is trait/risk increasing)',
#     'OR': 'Odds ratio (1 --> no effect; above 1 --> A1 is risk increasing)',
#     'BETA': '[linear/logistic] regression coefficient (0 --> no effect; above 0 --> A1 is trait/risk increasing)',
#     'LOG_ODDS': 'Log odds ratio (0 --> no effect; above 0 --> A1 is risk increasing)',
#     'INFO': 'INFO score (imputation quality; higher --> better imputation)',
#     'FRQ': 'Allele frequency',
#     'SIGNED_SUMSTAT': 'Directional summary statistic as specified by --signed-sumstats.',
#     'NSTUDY': 'Number of studies in which the SNP was genotyped.'
# }
#
# numeric_cols = ['P', 'N', 'N_CAS', 'N_CON', 'Z', 'OR', 'BETA', 'LOG_ODDS', 'INFO', 'FRQ', 'SIGNED_SUMSTAT', 'NSTUDY']
# "))

colmapping = c(
  SNP = "SNP",
  MARKERNAME = "SNP",
  SNPID = "SNP",
  RS = "SNP",
  RSID = "SNP",
  RS_NUMBER = "SNP",
  RS_NUMBERS = "SNP",
  NSTUDY = "NSTUDY",
  N_STUDY = "NSTUDY",
  NSTUDIES = "NSTUDY",
  N_STUDIES = "NSTUDY",
  P = "P",
  PVALUE = "P",
  P_VALUE = "P",
  PVAL = "P",
  P_VAL = "P",
  `P-VALUE` = "P",
  GC_PVALUE = "P",
  A1 = "A1",
  ALLELE1 = "A1",
  ALLELE_1 = "A1",
  EFFECT_ALLELE = "A1",
  REFERENCE_ALLELE = "A1",
  INC_ALLELE = "A1",
  EA = "A1",
  A2 = "A2",
  ALLELE2 = "A2",
  ALLELE_2 = "A2",
  OTHER_ALLELE = "A2",
  NON_EFFECT_ALLELE = "A2",
  DEC_ALLELE = "A2",
  NEA = "A2",
  N = "N",
  NCASE = "N_CAS",
  CASES_N = "N_CAS",
  N_CASE = "N_CAS",
  N_CASES = "N_CAS",
  N_CONTROLS = "N_CON",
  N_CAS = "N_CAS",
  N_CON = "N_CON",
  NCONTROL = "N_CON",
  CONTROLS_N = "N_CON",
  N_CONTROL = "N_CON",
  WEIGHT = "N",
  ZSCORE = "Z",
  `Z-SCORE` = "Z",
  GC_ZSCORE = "Z",
  Z = "Z",
  OR = "OR",
  B = "BETA",
  BETA = "BETA",
  STDERR = "SE",
  LOG_ODDS = "LOG_ODDS", EFFECTS = "BETA", EFFECT = "BETA",
  SIGNED_SUMSTAT = "SIGNED_SUMSTAT",
  INFO = "INFO", EAF = "FRQ", FRQ = "FRQ", MAF = "FRQ", FRQ_U = "FRQ",
  F_U = "FRQ",POS = "BP", BP = "BP")

numeric_cols = c("P", "N", "N_CAS", "N_CON", "Z", "OR", "BETA", "LOG_ODDS",
  "INFO", "FRQ", "SIGNED_SUMSTAT", "NSTUDY")

Reformat = function(sumstat){

sumstat = sumstat %>%
  rename_with(~clean_header(.x)) %>%
  rename_with(.cols = one_of(names(colmapping)),
              .fn = function(x){colmapping[x]})

sumstat = sumstat %>% mutate_if(
  names(.) %in% numeric_cols , as.numeric)

sumstat = sumstat %>%
  mutate(A1 = toupper(A1), A2 = toupper(A2))

return(sumstat)
}

complement <- function(x){
  switch (x,
          "A" = "T",
          "C" = "G",
          "T" = "A",
          "G" = "C",
          return(NA)
  )
}
complement = setNames(c("A", "C", "G", "T"), c("T", "G", "C", "A"))

#complement = c("A" = "T", "C" = "G","G" = "C", "T" = "A")

Align_allele =  function(sumstat.merge){

  flip.snps <- sumstat.merge %>% dplyr::filter(
    complement[A1.x] == A1.y, complement[A2.x] == A2.y) %>%
    dplyr::select(SNP)

  palindromic.snps <- sumstat.merge %>% dplyr::filter(
    (A1.x %in% c("A","T") &  A1.y %in% c("A","T")) |
      (A1.x %in% c("C","G") &  A1.y %in% c("C","G"))  ) %>%
    dplyr::select(SNP)

  recode.snps = sumstat.merge %>% dplyr::filter(
    A1.x == A2.y, A1.y == A2.x) %>%
    dplyr::select(SNP)

  const.snps = sumstat.merge %>% dplyr::filter(
    A1.x == A1.y, A2.x == A2.y) %>%
    dplyr::select(SNP)

  recode.flip.snps = sumstat.merge %>% dplyr::filter(
    complement[A1.x] == A2.y, complement[A2.x] == A1.y) %>%
    dplyr::select(SNP)

  return(list( flip.snps = setdiff(flip.snps$SNP, palindromic.snps$SNP),
               palindromic.snps = palindromic.snps$SNP,
               recode.snps = setdiff( recode.snps$SNP, palindromic.snps$SNP),
               recode.flip.snps = setdiff(recode.flip.snps$SNP, palindromic.snps$SNP),
               const.snps = const.snps$SNP ))
}


getInput <- function(sumstat1,
                     sumstat2,
                     marker.p.thres = 5e-8){


  # sumstat thresholding
  sumstat1 = sumstat1 %>% dplyr::filter(P < marker.p.thres)
  sumstat2 = sumstat2 %>% dplyr::filter(P < marker.p.thres)
  # merge snp sumstat dta
  sumstat.merge = inner_join(sumstat1, sumstat2, by = "SNP", suffix = c(".x", ".y"))
  snps.catag =  Align_allele(sumstat.merge)
  SNPs.kept = c(snps.catag$flip.snps, snps.catag$const.snps,
                snps.catag$recode.snps, snps.catag$recode.flip.snps)
  sumstat.merge = sumstat.merge %>% dplyr::filter(SNP %in% SNPs.kept)
  flip.sign.snp = c(snps.catag$recode.snps, snps.catag$recode.flip.snps)
  flip.sign = rep(1, nrow(sumstat.merge))
  flip.sign[sumstat.merge$SNP %in% flip.sign.snp] = -1
  sumstat.merge$BETA.y = flip.sign *  sumstat.merge$BETA.y

  sumstat.merge = sumstat.merge %>% dplyr::select(-c("A1.y","A2.y"))

  return(sumstat.merge)


}



# rename the summary statistics set

# change the minor allele frequency


r.to.z <- function(r, n){
  return(r * sqrt(n-2) /sqrt(1- r^2))
}

