# extract trait.names and N for the munge function
Aptamers.path <- 'D:\\Aptamers'
files <- list.files(Aptamers.path, pattern='.txt.csv.last')
aptamer.trait.names <- NULL
aptamers.N <- NULL
for (f in seq_along(files)) {
  aptamer.trait.names <- c(aptamer.trait.names, strsplit(files[f], '\\.')[[1]][1])
  LINE <- readLines(con = paste0(Aptamers.path, "\\", files[f]), n = 2)[-1]
  aptamers.N <- c(aptamers.N, as.numeric(tail(strsplit(LINE, '\t')[[1]], 1)))
}

# munge the summary statistics
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

# ldsc
Aptamers.munged.path <- '/Volumes/Elements/Aptamers_sumstat'
Aptamers.munged.files <- list.files(path=Aptamers.munged.path, pattern='.sumstats.gz', full.names=TRUE)[1:10]
Aptamer.trait.names <- NULL
files <- list.files(path=Aptamers.munged.path, pattern='.sumstats.gz')[1:10]
for(f in seq_along(files)) {Aptamer.trait.names <- c(Aptamer.trait.names, strsplit(files[f], '\\.')[[1]][1])}
aptamers.ldsc <- ldsc_parallel(
  traits = Aptamers.munged.files,
  ld = "~/R-workspace/eur_w_ld_chr",
  wld = "~/R-workspace/eur_w_ld_chr",
  sample.prev = rep(NA, length(Aptamer.trait.names)),
  population.prev = rep(NA, length(Aptamer.trait.names)),
  trait.names = paste0("aptamer_", Aptamer.trait.names),
  sep_weights = FALSE,
  chr = 22,
  n.blocks = 200,
  ldsc.log = 'aptamers',
  stand = TRUE,
  select = FALSE,
  chisq.max = NA,
  save.output = TRUE,
  output.path = 'data/',
  cores = 2
)
aptamers.ldsc.snowfall <- ldsc_snowfall(
  traits = Aptamers.munged.files,
  ld = "~/R-workspace/eur_w_ld_chr",
  wld = "~/R-workspace/eur_w_ld_chr",
  sample.prev = rep(NA, length(Aptamer.trait.names)),
  population.prev = rep(NA, length(Aptamer.trait.names)),
  trait.names = paste0("aptamer_", Aptamer.trait.names),
  sep_weights = FALSE,
  chr = 22,
  n.blocks = 200,
  ldsc.log = 'aptamers',
  stand = TRUE,
  select = FALSE,
  chisq.max = NA,
  save.output = TRUE,
  output.path = 'data/',
  cores = 2
)

require(Matrix)
rownames(aptamers.ldsc$S) <- colnames(aptamers.ldsc$S)
aptamers.ldsc$Ssmooth <- as.matrix((nearPD(x=aptamers.ldsc$S, corr = FALSE))$mat)

# Smooth the anthro[["S_Stand"]] matrix for EFA using the nearPD function in the Matrix package
aptamers.ldsc[["S_Standsmooth"]] <- as.matrix((nearPD(aptamers.ldsc[["S_Stand"]], corr = TRUE))$mat)

require(corpcor)
# Judge whether anthro[["S"]] is positive definite or not
is.positive.definite(aptamers.ldsc[["S"]])
# Partial Correlations Analysis
aptamers.ldsc[["pcorMatirx"]] <- cor2pcor(aptamers.ldsc[["S_Standsmooth"]])

require(psych)
# KMO test
KMO(r=aptamers.ldsc[["S_Standsmooth"]])
# Bartlett test
cortest.bartlett(R=aptamers.ldsc[["S_Standsmooth"]])

# Run the`paLDSC` function
paLDSC(S=aptamers.ldsc$S_Stand, V=aptamers.ldsc$V_Stand, r=100)

# Exploratory Factor Analysis (EFA) using ML/PA
# Run EFA with promax rotation and 2 factors using the factanal function in the stats package
aptamers.ldsc[["Ssmooth"]] <- as.matrix((nearPD(aptamers.ldsc[["S"]], corr = FALSE))$mat)
EFA.ml <- factanal(covmat = aptamers.ldsc[["Ssmooth"]], factors=2, rotation="promax")
EFA.mle <- fa(r=aptamers.ldsc[["Ssmooth"]], fm='mle', nfactors=2, rotate="promax", smooth=TRUE)
# In the psych package, run EFA.pa with promax rotation and 2-factors using the fa function.
# By default, improper correlations matrices are smoothed. This is not necessary for factor extraction using fm="pa" or fm="minres", but is needed to give some of the goodness of fit tests.
EFA.pa <- fa(r=aptamers.ldsc[["S"]], fm='pa', nfactors=2, rotate="promax", smooth=FALSE)
EFA.minres <- fa(r=aptamers.ldsc[["S"]], fm='minres', nfactors=2, rotate="promax", smooth=FALSE)

# Exploratory Factor Analysis (EFA) using spectral clustering
require(kernlab)
kernelMatrix <- as.kernelMatrix(as.matrix(aptamers.ldsc$Ssmooth), center=FALSE)
sc <- specc(x=kernelMatrix, centers=6, data=NULL, na.action = na.omit)

# Confirm the number of latent variables using hierarchical clustering
# In stats package, hclust method optional parameter: ward.D, ward.D2, single, complete, average, mcquitty, median, centroid
require(stats)
hclustTree <- hclust(d = as.dist((1 - aptamers.ldsc$S) / 2), method = "complete")
plot(hclustTree,
     main = "Hierarchical clustering for 87 GWAS of plasma protein levels measured with aptamers",
     xlab = 'Aptamers',
     ylab = 'Similarity [i.o. (1-Cor_coef)/2]',
     hang = -1)

# In pheatmap package, pheatmap clustering_method parameter is similar to the above stats::hclust method parameter option.
require(pheatmap)
pheatmap(
  aptamers.ldsc$S,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  clustering_method = "complete",
  annotation_row = data.frame(Aptamers = paste0("F", sc@.Data), row.names = rownames(aptamers.ldsc$S)),
  annotation_col = data.frame(Aptamers = paste0("F", sc@.Data), row.names = colnames(aptamers.ldsc$S)),
  main = "Hierarchical clustering and spectral clustering for genetic correlation matrix \n of 87 GWAS of plasma protein levels measured with aptamers",
  fontsize = 8,
  fontsize_row = 6,
  fontsize_col = 6
)

# Hierarchical Exploratory Factor Analysis (H-EFA) using BIRCH clustering
aptamers.ldsc$CFTree <- BirchCF(x=as.data.frame(aptamers.ldsc$S),
                                Type = 'df',
                                branchingfactor = 6,
                                threshold = 0.01)

# Specify the genomic confirmatory factor model
aptamers.ldsc$CFAofEFA <- write.model(S_LD=aptamers.ldsc$S, clusters=list(size=sc@size, .Data=sc@.Data), common = FALSE, hierarchical = FALSE, fix_resid = TRUE)

# Confirmatory factor analysis (CFA) based on specified the genomic confirmatory factor model
aptamers.ldsc$CFA <- usermodel(aptamers.ldsc, estimation = "DWLS", model = aptamers.ldsc$CFAofEFA, CFIcalc = TRUE, std.lv = TRUE, imp_cov = TRUE, fix_resid = TRUE, toler = NULL)
