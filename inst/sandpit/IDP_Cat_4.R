IDP.path <- '~/R-workspace/IDP/'
IDP.N <- c()

#munge the summary statistics
munge(
  files = paste0(IDP.path, '/', c(1326:1366), ".ma"),
  hm3 = "~/R-workspace/eur_w_ld_chr/w_hm3.snplist",
  trait.names = c(1326:1366),
  N = IDP.N,
  info.filter = 0.9,
  maf.filter = 0.01,
  column.names = list(SNP='SNP', A1='A1', A2='A2', effect='B', P='P', MAF='MAF'),
  parallel = TRUE,
  cores = 4,
  overwrite = TRUE
)

# ldsc
IDP.Cat.4 <- ldsc(
  traits = paste0(IDP.path, c(1326:1366), '.sumstats.gz'),
  ld = "~/R-workspace/eur_w_ld_chr",
  wld = "~/R-workspace/eur_w_ld_chr",
  sample.prev = rep(NA, length(c(1326:1366))),
  population.prev = rep(NA, length(c(1326:1366))),
  trait.names = c(1326:1366), 
  sep_weights = FALSE,
  chr = 22,
  n.blocks = 200,
  ldsc.log = 'IDP_Cat_4',
  stand = TRUE,
  select = FALSE,
  chisq.max = NA
)

# Convert genetic covariance matrix (i.e. S) to correlation matrix
require(Matrix)
rownames(IDP.Cat.4$S) <- colnames(IDP.Cat.4$S)
IDP.Cat.4$Ssmooth <- as.matrix((nearPD(x=IDP.Cat.4$S, corr = FALSE))$mat)
IDP.Cat.4.corMatirx <- as.matrix((nearPD(x=IDP.Cat.4$S, corr = TRUE))$mat)

# Exploratory Factor Analysis (EFA) using PCA and factor axis rotation
require(stats)
# Run EFA with promax rotation and 5 factors using the factanal function in the stats package
EFA <- factanal(covmat = IDP.Cat.4$Ssmooth, factors = 4, rotation = "promax")

# Confirm the number of latent variables using hierarchical clustering
# In stats package, hclust method optional parameter: ward.D, ward.D2, single, complete, average, mcquitty, median, centroid
require(stats)
hclustTree <- hclust(d = as.dist((1 - IDP.Cat.4.corMatirx) / 2), method = "complete")
plot(hclustTree,
     main = "Hierarchical clustering for 41 imaging-derived phenotypes (regional and tissue intensity)",
     xlab = 'Phenotypes',
     ylab = 'Similarity [i.o. (1-Cor_coef)/2]',
     hang = -1)

# Exploratory Factor Analysis (EFA) using spectral clustering
require(kernlab)
kernelMatrix <- as.kernelMatrix(as.matrix(IDP.Cat.4.corMatirx), center=FALSE)
sc <- specc(x=kernelMatrix, centers=5, data=NULL, na.action = na.omit)

# In pheatmap package, pheatmap clustering_method parameter is similar to the above stats::hclust method parameter option.
require(pheatmap)
pheatmap(
  IDP.Cat.4.corMatirx,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  clustering_method = "complete",
  annotation_row = data.frame(Traits = c(), row.names = rownames(IDP.Cat.4.corMatirx)),
  annotation_col = data.frame(Traits = c(), row.names = colnames(IDP.Cat.4.corMatirx)),
  main = "Hierarchical clustering and spectral clustering for genetic correlation matrix \n of 41 imaging-derived phenotypes (regional and tissue intensity)",
  fontsize = 8,
  fontsize_row = 6,
  fontsize_col = 6
)

# Hierarchical Exploratory Factor Analysis (H-EFA) using BIRCH clustering
IDP.Cat.4.CFTree <- BirchCF(x=as.data.frame(IDP.Cat.4.corMatirx),
                              Type = 'df',
                              branchingfactor = 6,
                              threshold = 0.01)

# Specify the genomic confirmatory factor model
CFAofEFA <- write.model(clusters=sc, S_LD=IDP.Cat.4$S, fix_resid=TRUE, common = FALSE, hierarchical = FALSE)

# Confirmatory factor analysis (CFA) based on specified the genomic confirmatory factor model
IDP.Cat.4.CFA <- usermodel(IDP.Cat.4, estimation = "DWLS", model = CFAofEFA, CFIcalc = TRUE, std.lv = TRUE, imp_cov = TRUE)
