#munge the summary statistics
munge(c("~/mrSEM/data/sumstats/GIANT/SNP_gwas_mc_merge_nogc.tbl.uniq",
        "~/mrSEM/data/sumstats/GIANT/GIANT_2015_WHR_COMBINED_EUR.txt",
        "~/mrSEM/data/sumstats/EGG/EGG_Obesity_Meta_Analysis_1.txt",
        "~/mrSEM/data/sumstats/GIANT/GIANT_2015_HIP_COMBINED_EUR.txt",
        "~/mrSEM/data/sumstats/GIANT/GIANT_2015_WC_COMBINED_EUR.txt",
        "~/mrSEM/data/sumstats/GIANT/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt",
        "~/mrSEM/data/sumstats/EGG/EGG_HC_DISCOVERY.v2.txt",
        "~/mrSEM/data/sumstats/EGG/EGG-GWAS-BL.txt",
        "~/mrSEM/data/sumstats/EGG/EGG_BW2_DISCOVERY.txt"),
      hm3 = "~/mrSEM/data/eur_w_ld_chr/w_hm3.snplist",
      trait.names = c("BMI2015", "waisthip", "childobese", "hip", "waist", "height", "headcirc", "birthlength", "birthweight"),
      N = c(NA, NA, 13848, NA, NA, NA, NA, NA, 26836),
      info.filter = 0.9,
      maf.filter = 0.01,
      # column.names = list(SNP=c('RSID', 'MarkerName'), A1='ref', A2='alt', effect='beta', P='pval', MAF='minor_AF', Z='tstat'),
      parallel = TRUE,
      cores = 8,
      overwrite = TRUE)

#run multivariate LDSC to create the S and V matrices
anthro <- ldsc(
  traits <- c(
    "~/mrSEM/data/sumstats/Anthro/BMI2015.sumstats.gz",
    "~/mrSEM/data/sumstats/Anthro/waisthip.sumstats.gz",
    "~/mrSEM/data/sumstats/Anthro/childobese.sumstats.gz",
    "~/mrSEM/data/sumstats/Anthro/waist.sumstats.gz",
    "~/mrSEM/data/sumstats/Anthro/hip.sumstats.gz",
    "~/mrSEM/data/sumstats/Anthro/height.sumstats.gz",
    "~/mrSEM/data/sumstats/Anthro/headcirc.sumstats.gz",
    "~/mrSEM/data/sumstats/Anthro/birthlength.sumstats.gz",
    "~/mrSEM/data/sumstats/Anthro/birthweight.sumstats.gz"
  ),
  sample.prev <- c(NA, NA, NA, NA, NA, NA, NA, NA, NA),
  population.prev <- c(NA, NA, NA, NA, NA, NA, NA, NA, NA),
  ld <- "~/mrSEM/data/eur_w_ld_chr/",
  wld <- "~/mrSEM/data/eur_w_ld_chr/",
  trait.names <- c("BMI", "WHR", "CO", "Waist", "Hip", "Height", "IHC", "BL", "BW")
)

require(Matirx)
# Convert genetic covariance matrix (i.e. S) to correlation matrix
anthro.corMatirx <- cov2cor(anthro[["S"]])
rownames(anthro.corMatirx) <- colnames(anthro[["S"]])

# Confirm the number of latent variables using hierarchical clustering
# In stats package, hclust method optional parameter: ward.D, ward.D2, single, complete, average, mcquitty, median, centroid
require(stats)
hclustTree <- stats::hclust(d = as.dist((1 - anthro.corMatirx) / 2), method = "complete")
plot(hclustTree,
     main = "Hierarchical clustering for 9 anthropometric traits",
     xlab = 'Traits',
     ylab = 'Similarity [i.o. (1-Cor_coef)/2]',
     hang = -1)

# Confirm the number of latent variables using hierarchical clustering
# In pheatmap package, pheatmap clustering_method parameter is similar to the above stats::hclust method parameter option.
require(pheatmap)
pheatmap::pheatmap(
  anthro.corMatirx,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  clustering_method = "complete",
  annotation_row = data.frame(Traits = c(rep("Overweight", 5), rep("Early life", 4)), row.names = rownames(anthro.corMatirx)),
  annotation_col = data.frame(Traits = c(rep("Overweight", 5), rep("Early life", 4)), row.names = colnames(anthro.corMatirx)),
  main = "Hierarchical clustering and spectral clustering for \n genetic correlation matrix of 9 anthropometric traits",
  fontsize = 8,
  fontsize_row = 6,
  fontsize_col = 6
)

# Exploratory Factor Analysis (EFA) using PCA and factor axis rotation
require(stats)
require(Matrix)
#In the Matrix package, smooth the anthro.corMatirx matrix for EFA using the nearPD function. 
Ssmooth <- as.matrix((nearPD(anthro.corMatirx, corr = FALSE))$mat)
#In the stats package, run EFA with promax rotation and 2 factors using the factanal function
EFA <- stats::factanal(covmat = Ssmooth, factors = 2, rotation = "promax")

# Exploratory Factor Analysis (EFA) using spectral clustering
require(kernlab)
kernelMatrix <- as.kernelMatrix(anthro.corMatirx, center=FALSE)
sc <- specc(x=kernelMatrix, centers=2, data=NULL, na.action = na.omit)

# Hierarchical Exploratory Factor Analysis (H-EFA) using BIRCH clustering
anthro.CFTree <-
  BirchCF(
    x = as.data.frame(anthro.corMatirx),
    Type = 'df',
    branchingfactor = 4,
    threshold = 0.4
  )

# Specify the genomic confirmatory factor model
CFAofEFA <- 'F1 =~ NA*BMI + WHR + CO + Waist + Hip
             F2 =~ NA*Height + IHC + BL + BW
             
             F1 ~~ F2
             Waist ~~ a*Waist
             a > .001'

# Confirmatory factor analysis (CFA) based on specified the genomic confirmatory factor model
anthro.CFA <- usermodel(anthro, estimation = "DWLS", model = CFAofEFA, CFIcalc = TRUE, std.lv = TRUE, imp_cov = TRUE)
