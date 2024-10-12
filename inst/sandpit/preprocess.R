# load variants.tsv
require(data.table)
variants <-
  fread(
    '~/mrSEM/data/sumstats/UKBB_GWAS_Imputed_v3/variants.tsv',
    select = c('variant', 'rsid', 'ref', 'alt'),
    data.table = FALSE
  )

# neuroticism.R: mood.txt, misery.txt, irrit.txt, hurt.txt, fedup.txt, nervous.txt, worry.txt, tense.txt, embarras.txt, nerves.txt, lonely.txt, guilt.txt
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

# ukb_pfactor.R
# Voices.txt
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
# Visions.txt
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
# Visions.txt
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
# Signs.txt
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
# Avoid.txt
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
# CannotRelax.txt
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
# Conspiracy.txt
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
# Depression.txt
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
# Foreboding.txt
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
# Inadequacy.txt
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
# Irritable.txt
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
# LostInterest.txt
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
# MultipleWorries.txt
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
# Mania.txt
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
# RecentWorry.txt
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
# RepeatedThoughts.txt
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
# Tired.txt
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
# Upset.txt
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
