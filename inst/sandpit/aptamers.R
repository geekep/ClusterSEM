
#extract trait.names and N for the munge function
Aptamers.path <- 'D:\\Aptamers'
files <- list.files(Aptamers.path)
aptamer.trait.names <- NULL
aptamers.N <- NULL
for (f in seq_along(files)) {
  aptamer.trait.names <- c(aptamer.trait.names, strsplit(files[f], '\\.')[[1]][1])
  LINE <- readLines(con = paste0(Aptamers.path, "\\", files[f]), n = 2)[-1]
  aptamers.N <- c(aptamers.N, as.numeric(tail(strsplit(LINE, '\t')[[1]], 1)))
}

#munge the summary statistics
munge(
  files = paste0(Aptamers.path, "\\", files),
  hm3 = "D:\\eur_w_ld_chr\\w_hm3.snplist",
  trait.names = aptamer.trait.names,
  N = aptamers.N,
  info.filter = 0.9,
  maf.filter = 0.01,
  log.name = 'Aptamers',
  column.names = list(SNP='SNP', A1='A1', A2='A2', effect='beta', P='p', MAF='EAF'),
  parallel = TRUE,
  cores = 3,
  overwrite = TRUE
)

#ldsc
aptamers.ldsc <- ldsc(
  traits = paste0(Aptamers.path, '\\', aptamer.trait.names, '.sumstats.gz'),
  ld = "D:\\eur_w_ld_chr",
  wld = "D:\\eur_w_ld_chr",
  sample.prev = rep(NA, length(aptamer.trait.names)),
  population.prev = rep(NA, length(aptamer.trait.names)),
  trait.names = paste0("aptamer_", aptamer.trait.names),
  sep_weights = FALSE,
  chr = 22,
  n.blocks = 200,
  ldsc.log = 'aptamers',
  stand = TRUE,
  select = FALSE,
  chisq.max = NA
)

# Convert genetic covariance matrix (i.e. S) to correlation matrix
require(Matrix)
rownames(aptamers.ldsc$S) <- colnames(aptamers.ldsc$S)
aptamers.ldsc$Ssmooth <- as.matrix((nearPD(x=aptamers.ldsc$S, corr = FALSE))$mat)
aptamers.ldsc$corMatirx <- as.matrix((nearPD(x=aptamers.ldsc$S, corr = TRUE))$mat)

# Exploratory Factor Analysis (EFA) using PCA and factor axis rotation
require(stats)
# Run EFA with promax rotation and 35-factors using the factanal function in the stats package
EFA <- factanal(covmat = aptamers.ldsc$Ssmooth, factors = 35, rotation = "promax")

# Confirm the number of latent variables using hierarchical clustering
# In stats package, hclust method optional parameter: ward.D, ward.D2, single, complete, average, mcquitty, median, centroid
require(stats)
hclustTree <- hclust(d = as.dist((1 - aptamers.ldsc$corMatirx) / 2), method = "complete")
plot(hclustTree,
     main = "Hierarchical clustering for 87 GWAS of plasma protein levels measured with aptamers",
     xlab = 'Aptamers',
     ylab = 'Similarity [i.o. (1-Cor_coef)/2]',
     hang = -1)

# Exploratory Factor Analysis (EFA) using spectral clustering
require(kernlab)
kernelMatrix <- as.kernelMatrix(as.matrix(aptamers.ldsc$corMatirx), center=FALSE)
sc <- specc(x=kernelMatrix, centers=6, data=NULL, na.action = na.omit)

# In pheatmap package, pheatmap clustering_method parameter is similar to the above stats::hclust method parameter option.
require(pheatmap)
pheatmap(
  aptamers.ldsc$corMatirx,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  clustering_method = "complete",
  annotation_row = data.frame(Aptamers = paste0("F", sc@.Data), row.names = rownames(aptamers.ldsc$corMatirx)),
  annotation_col = data.frame(Aptamers = paste0("F", sc@.Data), row.names = colnames(aptamers.ldsc$corMatirx)),
  main = "Hierarchical clustering and spectral clustering for genetic correlation matrix \n of 87 GWAS of plasma protein levels measured with aptamers",
  fontsize = 8,
  fontsize_row = 6,
  fontsize_col = 6
)

# Hierarchical Exploratory Factor Analysis (H-EFA) using BIRCH clustering
aptamers.ldsc.CFTree <- BirchCF(x=as.data.frame(aptamers.ldsc$corMatirx),
                            Type = 'df',
                            branchingfactor = 6,
                            threshold = 0.01)

# Specify the genomic confirmatory factor model
aptamers.ldsc.CFAofEFA <- write.model(S_LD=aptamers.ldsc$S, clusters=list(size=sc@size, .Data=sc@.Data), common = FALSE, hierarchical = FALSE, fix_resid = TRUE)

# Confirmatory factor analysis (CFA) based on specified the genomic confirmatory factor model
aptamers.CFA <- usermodel(aptamers.ldsc, estimation = "DWLS", model = aptamers.ldsc.CFAofEFA, CFIcalc = TRUE, std.lv = TRUE, imp_cov = TRUE, fix_resid = TRUE, toler = NULL)
