#munge the summary statistics
munge(
  files = c(
    "~/R-workspace/UKBB_GWAS_Imputed_v3/mood.txt",
    "~/R-workspace/UKBB_GWAS_Imputed_v3/misery.txt",
    "~/R-workspace/UKBB_GWAS_Imputed_v3/irrit.txt",
    "~/R-workspace/UKBB_GWAS_Imputed_v3/hurt.txt",
    "~/R-workspace/UKBB_GWAS_Imputed_v3/fedup.txt",
    "~/R-workspace/UKBB_GWAS_Imputed_v3/nervous.txt",
    "~/R-workspace/UKBB_GWAS_Imputed_v3/worry.txt",
    "~/R-workspace/UKBB_GWAS_Imputed_v3/tense.txt",
    "~/R-workspace/UKBB_GWAS_Imputed_v3/embarras.txt",
    "~/R-workspace/UKBB_GWAS_Imputed_v3/nerves.txt",
    "~/R-workspace/UKBB_GWAS_Imputed_v3/lonely.txt",
    "~/R-workspace/UKBB_GWAS_Imputed_v3/guilt.txt"
  ),
  hm3 = "~/R-workspace/eur_w_ld_chr/w_hm3.snplist",
  trait.names = c(
    "mood",
    "misery",
    "irrit",
    "hurt",
    "fedup",
    "nervous",
    "worry",
    "tense",
    "embarass",
    "nerves",
    "lonely",
    "guilt"
  ),
  N = c(
    329428,
    331856,
    322668,
    327832,
    330549,
    328725,
    328717,
    327232,
    323766,
    325248,
    332263,
    328769
  ), 
  info.filter = 0.9,
  maf.filter = 0.01,
  column.names = list(SNP='rsid', A1='ref', A2='alt', effect='beta', P='pval', MAF='minor_AF', Z='tstat'),
  parallel = TRUE,
  cores = 8,
  overwrite = TRUE
)

# ldsc
neuroticism <- ldsc(
  traits =
    c(
      "~/R-workspace/UKBB_GWAS_Imputed_v3/mood.sumstats.gz",
      "~/R-workspace/UKBB_GWAS_Imputed_v3/misery.sumstats.gz",
      "~/R-workspace/UKBB_GWAS_Imputed_v3/irrit.sumstats.gz",
      "~/R-workspace/UKBB_GWAS_Imputed_v3/hurt.sumstats.gz",
      "~/R-workspace/UKBB_GWAS_Imputed_v3/fedup.sumstats.gz",
      "~/R-workspace/UKBB_GWAS_Imputed_v3/nervous.sumstats.gz",
      "~/R-workspace/UKBB_GWAS_Imputed_v3/worry.sumstats.gz",
      "~/R-workspace/UKBB_GWAS_Imputed_v3/tense.sumstats.gz",
      "~/R-workspace/UKBB_GWAS_Imputed_v3/embarass.sumstats.gz",
      "~/R-workspace/UKBB_GWAS_Imputed_v3/nerves.sumstats.gz",
      "~/R-workspace/UKBB_GWAS_Imputed_v3/lonely.sumstats.gz",
      "~/R-workspace/UKBB_GWAS_Imputed_v3/guilt.sumstats.gz"
    ),
  ld = "~/R-workspace/eur_w_ld_chr",
  wld = "~/R-workspace/eur_w_ld_chr",
  sample.prev = c(
    .451,
    .427,
    .28,
    .556,
    .406,
    .237,
    .568,
    .171,
    .478,
    .213,
    .177,
    .283
  ),
  population.prev = c(
    .451,
    .427,
    .28,
    .556,
    .406,
    .237,
    .568,
    .171,
    .478,
    .213,
    .177,
    .283
  ),
  trait.names =
    c(
      "mood",
      "misery",
      "irritability",
      "hurt",
      "fedup",
      "nervous",
      "worry",
      "tense",
      "embarass",
      "nerves",
      "lonely",
      "guilt"
    ),
  sep_weights = FALSE,
  chr = 22,
  n.blocks = 200,
  ldsc.log = 'neuroticism',
  stand = TRUE,
  select = FALSE,
  chisq.max = NA
)

require(Matirx)
# Convert genetic covariance matrix (i.e. S) to correlation matrix
rownames(neuroticism[["S"]]) <- colnames(neuroticism[["S"]])
neuroticism.corMatirx <- cov2cor(neuroticism[["S"]])

# Confirm the number of latent variables using hierarchical clustering
# In stats package, hclust method optional parameter: ward.D, ward.D2, single, complete, average, mcquitty, median, centroid
require(stats)
hclustTree <- hclust(d = as.dist((1 - neuroticism.corMatirx) / 2), method = "complete")
plot(hclustTree,
     main = "Hierarchical clustering for 12 neuroticism traits",
     xlab = 'Traits',
     ylab = 'Similarity [i.o. (1-Cor_coef)/2]',
     hang = -1)

# Exploratory Factor Analysis (EFA) using PCA and factor axis rotation
require(Matrix)
#In the Matrix package, smooth the neuroticism.corMatirx matrix for EFA using the nearPD function. 
Ssmooth <- as.matrix((nearPD(neuroticism$S, corr = FALSE))$mat)
require(stats)
#In the stats package, run EFA with promax rotation and 3 factors using the factanal function 
EFA <- factanal(covmat = Ssmooth, factors = 3, rotation = "promax")

# Exploratory Factor Analysis (EFA) using spectral clustering
require(kernlab)
kernelMatrix <- as.kernelMatrix(as.matrix(neuroticism.corMatirx), center=FALSE)
sc <- specc(x=kernelMatrix, centers=3, data=NULL, na.action = na.omit)

# In pheatmap package, pheatmap clustering_method parameter is similar to the above stats::hclust method parameter option.
require(pheatmap)
pheatmap(
  neuroticism.corMatirx,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  clustering_method = "complete",
  annotation_row = data.frame(Traits = c("F1", "F1", "F1", "F2", "F1", "F3", "F3", "F3", "F2", "F3", "F1", "F2"), row.names = rownames(neuroticism.corMatirx)),
  annotation_col = data.frame(Traits = c("F1", "F1", "F1", "F2", "F1", "F3", "F3", "F3", "F2", "F3", "F1", "F2"), row.names = colnames(neuroticism.corMatirx)),
  main = "Hierarchical clustering and spectral clustering for \n genetic correlation matrix of 12 neuroticism traits",
  fontsize = 8,
  fontsize_row = 6,
  fontsize_col = 6
)

# Hierarchical Exploratory Factor Analysis (H-EFA) using BIRCH clustering
neuroticism.CFTree <- BirchCF(x=as.data.frame(neuroticism$S),
                               Type = 'df',
                               branchingfactor = 6,
                               threshold = 0.01)

# Specify the genomic confirmatory factor model
CFAofEFA <- 'F1 =~ NA*mood + misery + fedup +lonely + irritability
             F2 =~ NA*hurt + guilt + embarass
             F3 =~ NA*tense + nerves + nervous + worry'

# Confirmatory factor analysis (CFA) based on specified the genomic confirmatory factor model
neuroticism.CFA <- usermodel(neuroticism, estimation = "DWLS", model = CFAofEFA, CFIcalc = TRUE, std.lv = TRUE, imp_cov = TRUE, fix_resid = TRUE, toler = NULL)
