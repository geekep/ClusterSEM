# Lets start with a miniature example of 4 related variables.
# the first (t1) and the second (t2) causes the third (t3),
# which in turn cases the forth (t4). This means that while all 4 variables correlate,
# we can figure out that t1 doesn’t directly influence t4, by conditioning on t3.
t1 <- rnorm(1000)
t2 <- rnorm(1000)
t3 <- .4*t1 + .5*t2  + rnorm(1000)
t4 <- .6*t3 + rnorm(1000)

dataset <- cbind.data.frame(t1,t2,t3,t4)
# all variables are dependent
cor <- cor(dataset)

# t4, conditional on t2 is independent of t1
reg <- lm(t4 ~ t3 + t1, data=dataset)
summary(reg)

n <- 1000
network <- pc(suffStat = list(C = cor, n = n),
              indepTest = gaussCItest,
              labels = c("t1", "t2", "t3", "t4"),
              alpha = 0.01 )
plot(network, main="Toy DAG")

dev.off()

# Toy GWAS DAG
munge(files = c("~/R-workspace/mrSEM/data/ToyDAG/one.sums",
                "~/R-workspace/mrSEM/data/ToyDAG/two.sums",
                "~/R-workspace/mrSEM/data/ToyDAG/three.sums",
                "~/R-workspace/mrSEM/data/ToyDAG/four.sums",
                "~/R-workspace/mrSEM/data/ToyDAG/five.sums",
                "~/R-workspace/mrSEM/data/ToyDAG/six.sums") ,
      hm3 = "~/R-workspace/mrSEM/data/eur_w_ld_chr/w_hm3.snplist",
      N = rep(40000, 6),
      trait.names = c("one", "two", "three", "four", "five", "six"),
      info.filter = .9,
      maf.filter = .05)

covstruct <- ldsc(traits = c("~/R-workspace/mrSEM/inst/extdata/one.sumstats.gz",
                             "~/R-workspace/mrSEM/inst/extdata/two.sumstats.gz",
                             "~/R-workspace/mrSEM/inst/extdata/three.sumstats.gz",
                             "~/R-workspace/mrSEM/inst/extdata/four.sumstats.gz",
                             "~/R-workspace/mrSEM/inst/extdata/five.sumstats.gz",
                             "~/R-workspace/mrSEM/inst/extdata/six.sumstats.gz"),
                  trait.names = c("one", "two", "three", "four", "five", "six"),
                  ld = "~/R-workspace/mrSEM/data/eur_w_ld_chr/",
                  wld = "~/R-workspace/mrSEM/data/eur_w_ld_chr/",
                  sample.prev = NA,
                  population.prev = NA)

# Combining GenomicSEM and Bayesian network learning
suffstats <- GSEM.suffStat(suffStat = covstruct, rep = 1000)
pc_network <- pc(suffstats, indepTest = GSEM.indepTest, alpha = .05,
                 labels = c("one","two","three","four","five","six"))
plot(pc_network, main="simulated genomicSEM DAG")

