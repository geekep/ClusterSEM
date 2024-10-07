require(data.table)
require(dplyr)
require(stringr)
require(pheatmap)

# load variants.tsv
variants <-
  fread(
    '~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/variants.tsv',
    select = c('variant', 'rsid', 'ref', 'alt'),
    data.table = FALSE
  )
# Voices
Voices <-
  fread(
    '~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/Voices_20463.gwas.imputed_v3.both_sexes.tsv',
    data.table = FALSE
  )
Voices <- inner_join(Voices, variants, by = 'variant')
write.table(
  Voices,
  file = '~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/voices.txt',
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)
rm(Voices)
# Visions
Visions <-
  fread(
    '~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/Visions_20471.gwas.imputed_v3.both_sexes.tsv',
    data.table = FALSE
  )
Visions <- inner_join(Visions, variants, by = 'variant')
write.table(
  Visions,
  file = '~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/visions.txt',
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)
rm(Visions)
# Visions
Visions <-
  fread(
    '~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/Visions_20471.gwas.imputed_v3.both_sexes.tsv',
    data.table = FALSE
  )
Visions <- inner_join(Visions, variants, by = 'variant')
write.table(
  Visions,
  file = '~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/visions.txt',
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)
rm(Visions)
# Signs
Signs <-
  fread(
    '~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/Signs_20474.gwas.imputed_v3.both_sexes.tsv',
    data.table = FALSE
  )
Signs <- inner_join(Signs, variants, by = 'variant')
write.table(
  Signs,
  file = '~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/signs.txt',
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)
rm(Signs)
# Avoid
Avoid <-
  fread(
    '~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/Avoid_20495.gwas.imputed_v3.both_sexes.tsv',
    data.table = FALSE
  )
Avoid <- inner_join(Avoid, variants, by = 'variant')
write.table(
  Avoid,
  file = '~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/avoid.txt',
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)
rm(Avoid)
# CannotRelax
CannotRelax <-
  fread(
    '~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/CannotRelax_20515.gwas.imputed_v3.both_sexes.tsv',
    data.table = FALSE
  )
CannotRelax <- inner_join(CannotRelax, variants, by = 'variant')
write.table(
  CannotRelax,
  file = '~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/cannotRelax.txt',
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)
rm(CannotRelax)
# Conspiracy
Conspiracy <-
  fread(
    '~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/Conspiracy_20468.gwas.imputed_v3.both_sexes.tsv',
    data.table = FALSE
  )
Conspiracy <- inner_join(Conspiracy, variants, by = 'variant')
write.table(
  Conspiracy,
  file = '~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/conspiracy.txt',
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)
rm(Conspiracy)
# Depression
Depression <-
  fread(
    '~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/Depression.20510.gwas.imputed_v3.both_sexes.tsv',
    data.table = FALSE
  )
Depression <- inner_join(Depression, variants, by = 'variant')
write.table(
  Depression,
  file = '~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/depression.txt',
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)
rm(Depression)
# Foreboding
Foreboding <-
  fread(
    '~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/Foreboding_20512.gwas.imputed_v3.both_sexes.tsv',
    data.table = FALSE
  )
Foreboding <- inner_join(Foreboding, variants, by = 'variant')
write.table(
  Foreboding,
  file = '~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/foreboding.txt',
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)
rm(Foreboding)
# Inadequacy
Inadequacy <-
  fread(
    '~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/Inadequacy_20507.gwas.imputed_v3.both_sexes.tsv',
    data.table = FALSE
  )
Inadequacy <- inner_join(Inadequacy, variants, by = 'variant')
write.table(
  Inadequacy,
  file = '~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/inadequacy.txt',
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)
rm(Inadequacy)
# Irritable
Irritable <-
  fread(
    '~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/Irritable_20494.gwas.imputed_v3.both_sexes.tsv',
    data.table = FALSE
  )
Irritable <- inner_join(Irritable, variants, by = 'variant')
write.table(
  Irritable,
  file = '~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/irritable.txt',
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)
rm(Irritable)
# LostInterest
LostInterest <-
  fread(
    '~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/LostInterest_20514.gwas.imputed_v3.both_sexes.tsv',
    data.table = FALSE
  )
LostInterest <- inner_join(LostInterest, variants, by = 'variant')
write.table(
  LostInterest,
  file = '~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/lostInterest.txt',
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)
rm(LostInterest)
# MultipleWorries
MultipleWorries <-
  fread(
    '~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/MultipleWorries_20540.gwas.imputed_v3.both_sexes.tsv',
    data.table = FALSE
  )
MultipleWorries <- inner_join(MultipleWorries, variants, by = 'variant')
write.table(
  MultipleWorries,
  file = '~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/multipleWorries.txt',
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)
rm(MultipleWorries)
# Mania
Mania <-
  fread(
    '~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/Mania_20501.gwas.imputed_v3.both_sexes.tsv',
    data.table = FALSE
  )
Mania <- inner_join(Mania, variants, by = 'variant')
write.table(
  Mania,
  file = '~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/mania.txt',
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)
rm(Mania)
# RecentWorry
RecentWorry <-
  fread(
    '~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/RecentWorry_20520.gwas.imputed_v3.both_sexes.tsv',
    data.table = FALSE
  )
RecentWorry <- inner_join(RecentWorry, variants, by = 'variant')
write.table(
  RecentWorry,
  file = '~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/recentWorry.txt',
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)
rm(RecentWorry)
# RepeatedThoughts
RepeatedThoughts <-
  fread(
    '~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/RepeatedThoughts_20497.gwas.imputed_v3.both_sexes.tsv',
    data.table = FALSE
  )
RepeatedThoughts <- inner_join(RepeatedThoughts, variants, by = 'variant')
write.table(
  RepeatedThoughts,
  file = '~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/repeatedThoughts.txt',
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)
rm(RepeatedThoughts)
# Tired
Tired <-
  fread(
    '~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/Tired_20519.gwas.imputed_v3.both_sexes.tsv',
    data.table = FALSE
  )
Tired <- inner_join(Tired, variants, by = 'variant')
write.table(
  Tired,
  file = '~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/tired.txt',
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)
rm(Tired)
# Upset
Upset <-
  fread(
    '~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/Upset_20498.gwas.imputed_v3.both_sexes.tsv',
    data.table = FALSE
  )
Upset <- inner_join(Upset, variants, by = 'variant')
write.table(
  Upset,
  file = '~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/upset.txt',
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)
rm(Upset)

# munge
munge(
  files = c(
    "~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/voices.txt",
    "~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/visions.txt",
    "~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/signs.txt",
    "~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/conspiracy.txt",
    "~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/irritable.txt",
    "~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/mania.txt",
    "~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/depression.txt",
    "~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/inadequacy.txt",
    "~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/tired.txt",
    "~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/lostInterest.txt",
    "~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/avoid.txt",
    "~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/repeatedThoughts.txt",
    "~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/upset.txt",
    "~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/multipleWorries.txt",
    "~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/cannotRelax.txt",
    "~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/recentWorry.txt",
    "~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/foreboding.txt"
  ),
  hm3 = "~/mrSEM/data/eur_w_ld_chr/w_hm3.snplist",
  trait.names = c(
    "voices",
    "visions",
    "signs",
    "conspiracy",
    "irritable",
    "mania",
    "depression",
    "inadequacy",
    "tired",
    "lostInterest",
    "avoid",
    "repeatedThoughts",
    "upset",
    "multipleWorries",
    "cannotRelax",
    "recentWorry",
    "foreboding"
  ),
  N = NULL,
  info.filter = 0.9,
  maf.filter = 0.01,
  column.names = list(
    SNP = 'rsid',
    A1 = 'ref',
    A2 = 'alt',
    effect = 'beta',
    P = 'pval',
    MAF = 'minor_AF',
    Z = 'tstat',
    N = 'n_complete_samples'
  ), 
  parallel = TRUE,
  cores = 8,
  overwrite = TRUE
)

# ldsc
LDSCoutput_ukb_pfactor <- ldsc(
  traits =
    c(
      "~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/voices.sumstats.gz",
      "~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/visions.sumstats.gz",
      "~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/signs.sumstats.gz",
      "~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/conspiracy.sumstats.gz",
      "~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/irritable.sumstats.gz",
      "~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/mania.sumstats.gz",
      "~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/depression.sumstats.gz",
      "~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/inadequacy.sumstats.gz",
      "~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/tired.sumstats.gz",
      "~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/lostInterest.sumstats.gz",
      "~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/avoid.sumstats.gz",
      "~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/repeatedThoughts.sumstats.gz",
      "~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/upset.sumstats.gz",
      "~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/multipleWorries.sumstats.gz",
      "~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/cannotRelax.sumstats.gz",
      "~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/recentWorry.sumstats.gz",
      "~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/foreboding.sumstats.gz"
    ),
  ld = "~/mrSEM/data/eur_w_ld_chr",
  wld = "~/mrSEM/data/eur_w_ld_chr",
  sample.prev = rep(0.5, 17),
  population.prev = rep(0.5, 17),
  trait.names =
    c(
      "voices",
      "visions",
      "signs",
      "conspiracy",
      "irritable",
      "mania",
      "depression",
      "inadequacy",
      "tired",
      "lostInterest",
      "avoid",
      "repeatedThoughts",
      "upset",
      "multipleWorries",
      "cannotRelax",
      "recentWorry",
      "foreboding"
    ),
  sep_weights = FALSE,
  chr = 22,
  n.blocks = 200,
  ldsc.log = NULL,
  stand = FALSE,
  select = FALSE,
  chisq.max = NA
)

pheatmap(
  LDSCoutput_ukb_pfactor$S,
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete"
)

ukb_pfactor.CFTree <- BirchCF(x=as.data.frame(LDSCoutput_ukb_pfactor$S), Type = 'df', branchingfactor = 6, threshold = 0.01)
