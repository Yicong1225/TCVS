library(CVXR)
library(knockoff)
library(GUniFrac)
library(cluster)
library(dirmult)
install.packages("TCVS_1.0.tar.gz", repos = NULL)
library(TCVS)
data(P_60, package = "TCVS")
dim(P)
data(throat.tree, package = "GUniFrac")
data(throat.otu.tab, package = "GUniFrac")
data(DirMultOutput, package = "TCVS")

#### variable settings ####
n = 200
p = 60
setting = 1
q = 0.05
method = "BIC"
loop_start_index = 1
seed = 2023
type = "Dirmult" # type to generate compositional data
normalizeMethod = "Rowsum"
beta.slack.factor = 1
pseudocount = 0.5
maxlam = 0.1
minlam = 1e-7
nlam = 20

# generate simulated data
X <- get.all.OTU(n, p, trans = 0, type, normalizeMethod, pseudocount, seed = seed)
# clr transformation
Z <- get.all.OTU(n, p, trans = 2, type, normalizeMethod, pseudocount, seed = seed)
y <- get.all.y(Z, setting = setting, beta.slack.factor, seed = seed)

result.TCVS <- TCVS(
  X = X,
  Z = Z,
  y = y,
  P = P,
  method = method,
  maxlam,
  minlam,
  nlam,
  q = q,
  seed = seed
)
selected.OTU = which(result.TCVS$S.TCVS!=0)
selected.OTU # selected OTUs using TCVS

# -----------------------
# Result analysis
# -----------------------

if (setting != 2) {
  case <- c(1, 2, 4, 5, 8, 9, 10, 11, 12, 18, 19, 20)
}
if (setting == 2) {
  case <- c(1, 8, 16, 26, 29, 34, 41, 44, 47, 49, 51, 53)
}

S.TCVS = result.TCVS$S.TCVS
S.TCVS.true <- S.TCVS[, case]
S.TCVS.false <- S.TCVS[, -case]
suffixes <- c(".TCVS")
N = 1

if (N == 1) {
  for (suffix in suffixes) {
    assign(paste("TP", suffix, sep = ""), mean(get(paste("S", suffix, ".true", sep = ""))))
    assign(paste("FP", suffix, sep = ""), mean(get(paste("S", suffix, ".false", sep = ""))))
  }
} else {
  for (suffix in suffixes) {
    assign(paste("TP", suffix, sep = ""), apply(get(paste("S", suffix, ".true", sep = "")), 1, mean))
    assign(paste("FP", suffix, sep = ""), apply(get(paste("S", suffix, ".false", sep = "")), 1, mean))
  }
}

TPR <- numeric(length(suffixes))
FPR <- numeric(length(suffixes))
if (N == 1) {
  for (i in 1:length(suffixes)) {
    suffix <- suffixes[i]
    TPR[i] <- mean(get(paste("TP", suffix, sep = "")))
    FPR[i] <- mean(get(paste("FP", suffix, sep = "")))
  }
} else {
  for (i in 1:length(suffixes)) {
    suffix <- suffixes[i]
    TPR[i] <- mean(get(paste("TP", suffix, sep = ""))[1:N])
    FPR[i] <- mean(get(paste("FP", suffix, sep = ""))[1:N])
  }
}

names(TPR) <- paste("TPR", suffixes, sep = "")
names(FPR) <- paste("FPR", suffixes, sep = "")
TPR
FPR

