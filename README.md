
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SAT

<!-- badges: start -->

<!-- badges: end -->

The goal of SAT is to …

## Installation

You can install the released version of SAT from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("SAT")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("JianqiaoWang/SAT")
#> Using github PAT from envvar GITHUB_PAT
#> Downloading GitHub repo JianqiaoWang/SAT@HEAD
#>          checking for file 'C:\Users\wangjq\AppData\Local\Temp\RtmpgLNTvr\remotes26074f7279a\JianqiaoWang-SAT-2e7ee22/DESCRIPTION' ...  v  checking for file 'C:\Users\wangjq\AppData\Local\Temp\RtmpgLNTvr\remotes26074f7279a\JianqiaoWang-SAT-2e7ee22/DESCRIPTION'
#>       -  preparing 'SAT':
#>    checking DESCRIPTION meta-information ...     checking DESCRIPTION meta-information ...   v  checking DESCRIPTION meta-information
#>       -  checking for LF line-endings in source and make files and shell scripts
#>       -  checking for empty or unneeded directories
#>      NB: this package now depends on R (>=        NB: this package now depends on R (>= 3.5.0)
#>        WARNING: Added dependency on R >= 3.5.0 because serialized objects in  serialize/load version 3 cannot be read in older versions of R.  File(s) containing such objects: 'SAT/supp/geneTable.Final.rds'
#> -  building 'SAT_0.0.0.9000.tar.gz'
#>      
#> 
#> Installing package into 'C:/Users/wangjq/Documents/R/win-library/4.0'
#> (as 'lib' is unspecified)
```

## Example

This is a basic example which shows you how to solve a common problem:
Independent sequence

``` r
library(SAT)
## basic example code
p <- 1*10^5; ## total number of pairs
d1 = 100; d2 = 100; d3 = 1;
X <- c(rep(0,p-(d1 + d2 + d3)),rep(1,d1),rep(2,d2),rep(3,d3));
mu.1 = rep(0, p); mu.2 = rep(0, p);
mu.1[X==1|X==3] = 4; mu.2[X==2|X==3] = 5;
Z1 <- mu.1 + rnorm(p,0,1);
Z2 <- mu.2 + rnorm(p,0,1);
T1 = abs(Z1);T2 = abs(Z2);
P1 = 2* (1 - pnorm( T1)); P2 = 2* (1 - pnorm(T2));
MinMaxP.simu(P1, P2)
#> Ignoring the dependence
#> $M
#> [1] 4.235802e-06
#> 
#> $pv.sparse
#> [1] 0.0003344124
#> 
#> $pv.dense
#> [1] 0.000392061
```

For the typical genetic applications, we consider the LD between the
sequences. It requires the additional information of the external LD
matrix.

``` r
pvalue =  SAT::MinMaxP.geno(sumstat.comb$P.gwas, sumstat.comb$P.gtex, ref.geno = geno, ref.bed = "../geno/GTEx_WGS_838Indiv_Freeze_phased_MAF_0.01_chr1_22_EU",
                            output.dir = output.dir, var.name = sumstat.comb$variant_id, block.thresh = 0.99)
```
