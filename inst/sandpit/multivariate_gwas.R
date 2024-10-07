require(data.table)
require(dplyr)
require(stringr)
require(pheatmap)
require(GenomicSEM)
require(ctv)
require(kernlab)

# load variants.tsv
variants <-
  fread(
    '~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/variants.tsv',
    select = c('variant', 'rsid', 'ref', 'alt'),
    data.table = FALSE
  )

# prepare mood.txt, misery.txt, irrit.txt, hurt.txt, fedup.txt, nervous.txt, worry.txt, tense.txt, embarras.txt, nerves.txt, lonely.txt, guilt.txt
# mood.txt
Mood <-
  fread(
    '~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/Mood_1920.gwas.imputed_v3.both_sexes.tsv',
    data.table = FALSE
  )
Mood <- inner_join(Mood, variants, by = 'variant')
write.table(
  Mood,
  file = '~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/mood.txt',
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)
rm(Mood)

# misery.txt
Misery <-
  fread(
    '~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/Misery_1930.gwas.imputed_v3.both_sexes.tsv',
    data.table = FALSE
  )
Misery <- inner_join(Misery, variants, by = 'variant')
write.table(
  Misery,
  file = '~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/misery.txt',
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)
rm(Misery)

# irrit.txt
Irritability <-
  fread(
    '~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/Irritability_1940.gwas.imputed_v3.both_sexes.tsv',
    data.table = FALSE
  )
Irritability <- inner_join(Irritability, variants, by = 'variant')
write.table(
  Irritability,
  file = '~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/irrit.txt',
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)
rm(Irritability)

# hurt.txt
Hurt <-
  fread(
    '~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/Hurt_1950.gwas.imputed_v3.both_sexes.tsv',
    data.table = FALSE
  )
Hurt <- inner_join(Hurt, variants, by = 'variant')
write.table(
  Hurt,
  file = '~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/hurt.txt',
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)
rm(Hurt)

# fedup.txt
Fedup <-
  fread(
    '~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/Fedup_1960.gwas.imputed_v3.both_sexes.tsv',
    data.table = FALSE
  )
Fedup <- inner_join(Fedup, variants, by = 'variant')
write.table(
  Fedup,
  file = '~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/fedup.txt',
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)
rm(Fedup)

# nervous.txt
Nervous <-
  fread(
    '~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/Nervous_1970.gwas.imputed_v3.both_sexes.tsv',
    data.table = FALSE
  )
Nervous <- inner_join(Nervous, variants, by = 'variant')
write.table(
  Nervous,
  file = '~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/nervous.txt',
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)
rm(Nervous)

# worry.txt
Worry <-
  fread(
    '~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/Worry_1980.gwas.imputed_v3.both_sexes.tsv',
    data.table = FALSE
  )
Worry <- inner_join(Worry, variants, by = 'variant')
write.table(
  Worry,
  file = '~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/worry.txt',
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)
rm(Worry)

# tense.txt
Tense <-
  fread(
    '~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/Tense_1990.gwas.imputed_v3.both_sexes.tsv',
    data.table = FALSE
  )
Tense <- inner_join(Tense, variants, by = 'variant')
write.table(
  Tense,
  file = '~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/tense.txt',
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)
rm(Tense)

# embarras.txt
Embarras <-
  fread(
    '~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/Embarras_2000.gwas.imputed_v3.both_sexes.tsv',
    data.table = FALSE
  )
Embarras <- inner_join(Embarras, variants, by = 'variant')
write.table(
  Embarras,
  file = '~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/embarras.txt',
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)
rm(Embarras)

# nerves.txt
Nerves <-
  fread(
    '~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/Nerves_2010.gwas.imputed_v3.both_sexes.tsv',
    data.table = FALSE
  )
Nerves <- inner_join(Nerves, variants, by = 'variant')
write.table(
  Nerves,
  file = '~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/nerves.txt',
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)
rm(Nerves)

# lonely.txt
Lonely <-
  fread(
    '~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/Lonely_2020.gwas.imputed_v3.both_sexes.tsv',
    data.table = FALSE
  )
Lonely <- inner_join(Lonely, variants, by = 'variant')
write.table(
  Lonely,
  file = '~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/lonely.txt',
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)
rm(Lonely)

# guilt.txt
Guilt <-
  fread(
    '~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/Guilt_2030.gwas.imputed_v3.both_sexes.tsv',
    data.table = FALSE
  )
Guilt <- inner_join(Guilt, variants, by = 'variant')
write.table(
  Guilt,
  file = '~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/guilt.txt',
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)
rm(Guilt)


# munge
munge(
  files = c(
    "~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/mood.txt",
    "~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/misery.txt",
    "~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/irrit.txt",
    "~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/hurt.txt",
    "~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/fedup.txt",
    "~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/nervous.txt",
    "~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/worry.txt",
    "~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/tense.txt",
    "~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/embarras.txt",
    "~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/nerves.txt",
    "~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/lonely.txt",
    "~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/guilt.txt"
  ),
  hm3 = "~/mrSEM/data/eur_w_ld_chr/w_hm3.snplist",
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
LDSCoutput_multivariate <- ldsc(
    traits =
      c(
        "~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/mood.sumstats.gz",
        "~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/misery.sumstats.gz",
        "~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/irrit.sumstats.gz",
        "~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/hurt.sumstats.gz",
        "~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/fedup.sumstats.gz",
        "~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/nervous.sumstats.gz",
        "~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/worry.sumstats.gz",
        "~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/tense.sumstats.gz",
        "~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/embarass.sumstats.gz",
        "~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/nerves.sumstats.gz",
        "~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/lonely.sumstats.gz",
        "~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/guilt.sumstats.gz"),
    ld = "~/mrSEM/data/eur_w_ld_chr",
    wld = "~/mrSEM/data/eur_w_ld_chr",
    sample.prev = c(0.451,0.427,0.28,0.556,.406,.237,.568,.171,.478,.213,.177,.283),
    population.prev = c(0.451,0.427,0.28,0.556,.406,.237,.568,.171,.478,.213,.177,.283),
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
    ldsc.log = NULL,
    stand = FALSE,
    select = FALSE,
    chisq.max = NA
  )

neuroticism.corMatirx <- cov2cor(LDSCoutput_multivariate[["S"]])
rownames(neuroticism.corMatirx) <- colnames(LDSCoutput_multivariate[["S"]])

# method: ward.D, ward.D2, single, complete, average, mcquitty, median, centroid
hclustTree <- hclust(d = as.dist((1 - neuroticism.corMatirx) / 2), method = "complete")
plot(hclustTree,
     main = "Hierarchical clustering for 12 neuroticism traits",
     xlab = 'Traits',
     ylab = '(1-Correlation_coefficient)/2',
     hang = -1)

# pheatmap
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

# Spectral Clustering
kernelMatrix <- as.kernelMatrix(as.matrix(neuroticism.corMatirx), center=FALSE)
sc <- specc(x=kernelMatrix, centers=3, data=NULL, na.action = na.omit)

# multivariate.CFTree <- BirchCF(x=as.data.frame(LDSCoutput_multivariate$S),
#                                Type = 'df',
#                                branchingfactor = 6,
#                                threshold = 0.01)

# smooth the neuroticism matrix for EFA using the nearPD function in the Matrix package. 
Ssmooth <- as.matrix((nearPD(LDSCoutput_multivariate$S, corr = FALSE))$mat)

#run EFA with promax rotation and 3 factors using the factanal function in the stats package
EFA <- factanal(covmat = Ssmooth, factors = 3, rotation = "promax")

#Specify the Genomic confirmatory factor model
CFAofEFA <- 'F1 =~ NA*mood + misery + fedup +lonely + irritability
             F2 =~ NA*hurt + guilt + embarass
             F3 =~ NA*tense + nerves + nervous + worry'

#run the model
Neuroticism <- usermodel(anthro, estimation = "DWLS", model = CFAofEFA, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)


