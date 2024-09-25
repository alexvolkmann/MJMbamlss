# New simulation scenario

library(MFPCA)
library(MJMbamlss)
library(gmfamm)
i <- 1500

# Sample Sizes
n <- c(50, 250, 500)

# Numbers of longitudinal markers
k <- c(5, 10, 15)

# Correlation between markers
r <- c(0.9, 0.5, 0.2)

# Marker sparsity
nk <- c(5, 10, 20)
argvals <- seq(0, 1, by = 0.01)

# Censoring rate
# TO BE DISCUSSED

# Number of basis functions
nb <- 3
for (N in n) {
  for (K in k) {

  }
}

B <- rep(list(funData::eFun(argvals, nb, type = "FourierLin")), k[1])
# autcor <- c(1, 0.8, 0.75)
# crocor <- c(0.9, 0.8, 0.75)
# a <- c(autcor,
#        r[1] * rep(1/(2:k[1]), each = length(crocor)) *
#          rep(crocor, k[1] - 1))

fun <- Vectorize(function (x, rho) exp(1/(rho)^2*(1/(15*3) - x/(15*3))),
                 vectorize.args = "x")
plot(1:(5*nb), fun(1:(5*nb), rho = 0.9))
eigen(toeplitz(fun(1:(15*nb), rho = 0.5)))$values
COV <- toeplitz(fun(1:(k[1]*nb), rho = r[1]))
mfpca <- MFPCA_cov(cov = COV, basis_funs = B)
plot(mfpca$functions, obs = 1:3)
cumsum(mfpca$values) / sum(mfpca$values)

psi <- gmfamm::simMuFu(type = "split", argvals = rep(list(argvals), K), M = 5,
                       eFunType = "Fourier", eValType = "linear", N = N,
                       seed = i)
