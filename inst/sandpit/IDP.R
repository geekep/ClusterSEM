require(GenomicSEM)
require(pheatmap)

files.name <- list.files('/storageB/xuzhen/IDP')
files.name.vector <- c()
trait.name.vector <- c()
N.vector <- c()
for(file.name in files.name)
{
  # files.name.vector <- c(files.name.vector, paste0('/storageB/xuzhen/IDP/',file.name))
  trait.name.vector <- c(trait.name.vector, substr(strsplit(file.name, '/')[[1]][length(strsplit(file.name, '/')[[1]])], 1, 4))
  # df <- read.table(file, head=TRUE)
  # N.vector <- c(N.vector, df$N[1])
  # rm(df)
}

# munge
munge(
  files = files.name.vector,
  hm3 = "~/mrSEM/data/eur_w_ld_chr/w_hm3.snplist",
  trait.names = trait.name.vector,
  N = N.vector,
  info.filter = 0.9,
  maf.filter = 0.01,
  column.names = list(SNP='SNP', A1='A1', A2='A2', effect='B', P='P', MAF='MAF'),
  parallel = TRUE,
  cores = 24,
  overwrite = TRUE
)

# ldsc
LDSCoutput_IDP_1326_1366 <- ldsc(
  traits = paste0("/storageB/xuzhen/IDP_sumstat/", trait.name.vector[1326:1366], '.sumstats.gz'),
  ld = "~/mrSEM/data/eur_w_ld_chr",
  wld = "~/mrSEM/data/eur_w_ld_chr",
  sample.prev = rep(0.5, length(trait.name.vector[1326:1366])),
  population.prev = rep(0.5, length(trait.name.vector[1326:1366])),
  trait.names = trait.name.vector[1326:1366], 
  sep_weights = FALSE,
  chr = 22,
  n.blocks = 200,
  ldsc.log = NULL,
  stand = FALSE,
  select = FALSE,
  chisq.max = NA
)

pheatmap(
  LDSCoutput_IDP$S,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "single"
)

IDP.CFTree.2 <- BirchCF(x=as.data.frame(LDSCoutput_IDP_1326_1366$S), Type = 'df', branchingfactor = 4, threshold = 0.4)
