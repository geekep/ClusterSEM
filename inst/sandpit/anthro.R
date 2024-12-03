# munge the summary statistics
munge(c("~/mrSEM/data/sumstats/GIANT/SNP_gwas_mc_merge_nogc.tbl.uniq",
        "~/mrSEM/data/sumstats/GIANT/GIANT_2015_WHR_COMBINED_EUR.txt",
        "~/mrSEM/data/sumstats/EGG/EGG_Obesity_Meta_Analysis_1.txt",
        "~/mrSEM/data/sumstats/GIANT/GIANT_2015_HIP_COMBINED_EUR.txt",
        "~/mrSEM/data/sumstats/GIANT/GIANT_2015_WC_COMBINED_EUR.txt",
        "~/mrSEM/data/sumstats/GIANT/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt",
        "~/mrSEM/data/sumstats/EGG/EGG_HC_DISCOVERY.v2.txt",
        "~/mrSEM/data/sumstats/EGG/EGG-GWAS-BL.txt",
        "~/mrSEM/data/sumstats/EGG/EGG_BW2_DISCOVERY.txt"),
      hm3 = "~/R-workspace/eur_w_ld_chr/w_hm3.snplist",
      trait.names = c("BMI2015", "waisthip", "childobese", "hip", "waist", "height", "headcirc", "birthlength", "birthweight"),
      N = c(NA, NA, 13848, NA, NA, NA, NA, NA, 26836),
      info.filter = 0.9,
      maf.filter = 0.01,
      # column.names = list(SNP=c('RSID', 'MarkerName'), A1='ref', A2='alt', effect='beta', P='pval', MAF='minor_AF', Z='tstat'),
      parallel = TRUE,
      cores = 8,
      overwrite = TRUE)

# run multivariate LDSC to create the S and V matrices
anthro <- ldsc(
  traits <- c(
    "~/R-workspace/Anthro/BMI2015.sumstats.gz",
    "~/R-workspace/Anthro/waisthip.sumstats.gz",
    "~/R-workspace/Anthro/childobese.sumstats.gz",
    "~/R-workspace/Anthro/waist.sumstats.gz",
    "~/R-workspace/Anthro/hip.sumstats.gz",
    "~/R-workspace/Anthro/height.sumstats.gz",
    "~/R-workspace/Anthro/headcirc.sumstats.gz",
    "~/R-workspace/Anthro/birthlength.sumstats.gz",
    "~/R-workspace/Anthro/birthweight.sumstats.gz"
  ),
  sample.prev <- c(NA, NA, NA, NA, NA, NA, NA, NA, NA),
  population.prev <- c(NA, NA, NA, NA, NA, NA, NA, NA, NA),
  ld <- "~/R-workspace/eur_w_ld_chr/",
  wld <- "~/R-workspace/eur_w_ld_chr",
  trait.names <- c("BMI", "WHR", "CO", "Waist", "Hip", "Height", "IHC", "BL", "BW"),
  sep_weights = FALSE,
  chr = 22,
  n.blocks = 200,
  ldsc.log = 'anthro',
  stand = TRUE,
  select = FALSE,
  chisq.max = NA
)

require(Matirx)
rownames(anthro[["S"]]) <- colnames(anthro[["S"]])
rownames(anthro[["S_Stand"]]) <- colnames(anthro[["S_Stand"]])

# Smooth the anthro[["S_Stand"]] matrix for EFA using the nearPD function in the Matrix package
anthro[["S_Standsmooth"]] <- as.matrix((nearPD(anthro[["S_Stand"]], corr = TRUE))$mat)

require(corpcor)
# Judge whether anthro[["S"]] is positive definite or not
is.positive.definite(anthro[["S"]])
# Partial Correlations Analysis
anthro[["pcorMatirx"]] <- cor2pcor(anthro[["S_Standsmooth"]])

require(psych)
# KMO test
KMO(r=anthro[["S_Standsmooth"]])
# Bartlett test
cortest.bartlett(R=anthro[["S_Standsmooth"]])

# Run the`paLDSC` function
paLDSC(S=anthro$S_Stand, V=anthro$V_Stand, r=100)

# Exploratory Factor Analysis (EFA) using ML/PA
# Run EFA with promax rotation and 2 factors using the factanal function in the stats package
anthro[["Ssmooth"]] <- as.matrix((nearPD(anthro[["S"]], corr = FALSE))$mat)
EFA.ml <- factanal(covmat = anthro[["Ssmooth"]], factors=2, rotation="promax")
EFA.mle <- fa(r=anthro[["Ssmooth"]], fm='mle', nfactors=2, rotate="promax", smooth=FALSE)
# In the psych package, run EFA.pa with promax rotation and 2-factors using the fa function.
# By default, improper correlations matrices are smoothed. This is not necessary for factor extraction using fm="pa" or fm="minres", but is needed to give some of the goodness of fit tests.
EFA.pa <- fa(r=anthro[["S"]], fm='pa', nfactors=2, rotate="promax", smooth=FALSE)
EFA.minres <- fa(r=anthro[["S"]], fm='minres', nfactors=2, rotate="promax", smooth=FALSE)

# Exploratory Factor Analysis (EFA) using spectral clustering
require(kernlab)
kernelMatrix <- as.kernelMatrix(anthro[["S_Stand"]], center=FALSE)
sc <- specc(x=kernelMatrix, centers=2, data=NULL, na.action = na.omit)

# Confirm the number of latent variables using hierarchical clustering
# In stats package, hclust method optional parameter: ward.D, ward.D2, single, complete, average, mcquitty, median, centroid
require(stats)
hclustTree <- hclust(d = as.dist((1 - anthro[["S_Stand"]]) / 2), method = "complete")
plot(hclustTree,
     main = "Hierarchical clustering for 9 anthropometric traits",
     xlab = 'Traits',
     ylab = 'Similarity [i.o. (1-Cor_coef)/2]',
     hang = -1)

# In pheatmap package, pheatmap clustering_method parameter is similar to the above stats::hclust method parameter option.
require(pheatmap)
pheatmap(
  anthro[["S_Stand"]],
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  clustering_method = "complete",
  annotation_row = data.frame(Traits = c(rep("Overweight", 5), rep("Early life", 4)), row.names = rownames(anthro[["S_Stand"]])),
  annotation_col = data.frame(Traits = c(rep("Overweight", 5), rep("Early life", 4)), row.names = colnames(anthro[["S_Stand"]])),
  main = "Hierarchical clustering and spectral clustering for \n genetic correlation matrix of 9 anthropometric traits",
  fontsize = 8,
  fontsize_row = 6,
  fontsize_col = 6
)

# Hierarchical Exploratory Factor Analysis (H-EFA) using BIRCH clustering
anthro[["CFTree"]] <- BirchCF(x = as.data.frame(anthro[["S_Stand"]]),
                              Type = 'df', branchingfactor = 4, threshold = 0.4)

# Specify the genomic confirmatory factor model
anthro[["CFAofEFA"]] <- 'F1 =~ NA*BMI + WHR + CO + Waist + Hip
                         F2 =~ NA*Height + IHC + BL + BW
                         F1 ~~ F2
                         Waist ~~ a*Waist
                         a > .001'

# Confirmatory factor analysis (CFA) based on specified the genomic confirmatory factor model
anthro[["CFA"]] <- usermodel(anthro, estimation = "DWLS", model = anthro[["CFAofEFA"]], CFIcalc = TRUE, std.lv = TRUE, imp_cov = TRUE, fix_resid = TRUE, toler = NULL)
