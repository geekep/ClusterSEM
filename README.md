# ClusterSEM (Structural Equation Models based on Clustering)

## Installation

``` {#Installation .R}
install.packages("devtools")
devtools::install_github('geekep/ClusterSEM')
```

## Usage

### Step 1: munge the summary statistics

``` r
munge(
  files,
  hm3,
  trait.names,
  N,
  info.filter = 0.9,
  maf.filter = 0.01,
  column.names = list(SNP='SNP', A1='A1', A2='A2', effect='B', P='P', MAF='MAF'),
  parallel = TRUE,
  cores = 4,
  overwrite = TRUE
)
```

### Step 2: Estimate genetic covariance structure (i.o., LDSC, HDL, metaCCA)

``` R
ldsc_output <- ldsc(
  traits,
  ld = "eur_w_ld_chr",
  wld = "eur_w_ld_chr",
  sample.prev,
  population.prev,
  trait.names, 
  sep_weights = FALSE,
  chr = 22,
  n.blocks = 200,
  ldsc.log = “”,
  stand = TRUE,
  select = FALSE,
  chisq.max = NA
)
```

### Step 3: Specify the genomic confirmatory factor model using spectral clustering

``` R
CFAofEFA <- write.model(S_LD=ldsc_output$S)
```

### Step 4: Confirmatory factor analysis (CFA) based on specified the genomic confirmatory factor model

``` R
usermodel(ldsc_output, estimation = "DWLS", model = CFAofEFA, CFIcalc = TRUE, std.lv = TRUE, imp_cov = TRUE, fix_resid = TRUE, toler = NULL)
```
