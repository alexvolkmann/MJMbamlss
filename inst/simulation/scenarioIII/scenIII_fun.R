
# Run a Simulation Scenario -----------------------------------------------


# Set up R session --------------------------------------------------------



# Specify location
results_wd <- paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.",
                     "hu-berlin.de,share=volkmale/MJMbamlss/simulation")


# Always
library(survival)
library(JMbayes2)
library(bamlss)
library(MFPCA)
library(tidyverse)
library(parallel)
library(Rcpp)
library(Matrix)
library(sparseFLMM)
library(MJMbamlss)



setting <- "scen_III"



# Different Specifications ------------------------------------------------

# Specifications are ordered so that estimation should get successively more
# computationally expensive as more data are available


# Number of individuals
N <- c(100, 250, 500)

# Numbers of longitudinal markers
K <- c(5, 10, 15)

# Correlation between markers
R <- c(0.9, 0.5, 0.2)

# Marker sparsity
NK <- c(0.9, 0.75, 0.5)

# Censoring rate
CE <- c(2, 3, 4)



# Objects for all clusters ------------------------------------------------

# Functional observation grid
argvals <- seq(0, 1, by = 0.01)

# Function to calculate the elements of the Toeplitz-Matrix
toeplitzfun <- Vectorize(function (x, rho) {
  exp(1 / (rho)^2 * (1/(20*3) - x/(20*3)))
}, vectorize.args = "x")

# Number of Fourier basis functions
nb <- 3


# Calculate and save MFPCs
dir.create(file.path(results_wd, setting, "mfpcs"),
           showWarnings = FALSE)
for (k in K) {
  for (r in R) {

    # Fourier basis functions
    B <- rep(list(funData::eFun(argvals, nb, type = "Fourier")), k)

    # Toeplitz covariance for basis coefficients
    COV <- toeplitz(toeplitzfun(1:(k*nb), rho = r))

    # Calculate MFPCs
    mfpca <- MFPCA_cov(cov = COV, basis_funs = B)

    saveRDS(mfpca,
            file = file.path(
              results_wd, setting, "mfpcs", paste0(
                paste("mfpc", k, r, sep = "_"), ".rds"
              )))
  }
}



# Simulate one scenario ---------------------------------------------------


sim_scen <- function (n, k, r, nk, ce, i) {

  # Name for scenario
  scen <- paste("set", n, k, r, nk, ce, sep = "_")

  # Create the folder
  dir.create(file.path(results_wd, setting, scen),
             showWarnings = FALSE)

  # Only create the data set if it is the first time
  # -- -- -- -- --
  if (! file.exists(file.path(results_wd, setting, scen, "data",
                              paste0("d", i, ".rds")))) {

    dir.create(file.path(results_wd, setting, scen, "data"),
               showWarnings = FALSE)


    # Load the corresponding MFPCs
    m <- readRDS(file.path(
      results_wd, setting, "mfpcs", paste0(
        paste("mfpc", k, r, sep = "_"), ".rds"
      )))

    set.seed(i)

    # Simulate the data
    d_sim <- simMultiJM(nsub = n, times = argvals,
                        max_obs = 15, probmiss = nk, maxfac = ce,
                        nmark = k, long_assoc = "FPC", M = k * nb,
                        FPC_bases = m$functions, FPC_evals = m$values,
                        ncovar = 2,
                        lambda = function(t, x) {
                          1.65 * t^(0.65)
                        },
                        gamma = function(x) {
                          -7 + 0.3*x[, 3]
                        },
                        alpha = rep(list(
                          function(t, x) 1,
                          function(t, x) 0.775,
                          function(t, x) 0.550,
                          function(t, x) 0.325,
                          function(t, x) 0.1), k %/% 5),
                        mu = rep(list(function(t, x, r){
                          0 + 1*t + 0.3*x[, 3] + 0.3*t*x[, 3]
                        }), k),
                        sigma = function(t, x) {
                          log(0.06) + 0*t
                        },
                        tmax = NULL, seed = NULL,
                        full = TRUE, file = NULL)

    # Uncomment to save simulated data
    saveRDS(d_sim,
            file = file.path(results_wd, setting, scen, "data",
                             paste0("d", i, ".rds")))

  }
  # -- -- -- -- --

  # Load the data set
  d_sim <- readRDS(file.path(results_wd, setting, scen, "data",
                             paste0("d", i, ".rds")))

  d <- d_sim$data
  mfpcs <- d_sim$fpc_base

  # Use all true data generating FPCs
  nfpc <- nObs(d_sim$fpc_base)
  mfpca_est_list <- lapply(seq_len(nfpc), function (it) {
    list(functions = extractObs(mfpcs, it),
         values = mfpcs$values[it])
  })

  # Model formula
  form <- list(
    Surv2(survtime, event, obs = y) ~ -1 +
      s(survtime, k = 20, bs = "ps"),
    gamma ~ 1 + x3,
    as.formula(paste0(
      "mu ~ -1 + marker + obstime:marker + x3:marker +",
      " obstime:x3:marker +",
      paste0(lapply(seq_len(nfpc), function(x) {
        paste0("s(id, fpc.", x, ", bs = 'unc_pcre', xt = list('mfpc' = ",
               "mfpca_est_list[[", x, "]], 'scale' = FALSE))")
      }), collapse = " + "))),
    sigma ~ -1 + marker,
    alpha ~ -1 + marker
  )

  t_est <- system.time(
    b_est <- bamlss(form, family = mjm_bamlss,
                    data = d,
                    timevar = "obstime", maxit = 1500, n.iter = 5500,
                    burnin = 500, thin = 5, verbose = TRUE)
  )
  attr(b_est, "comp_time") <- t_est
  attr(b_est, "FPCs") <- mfpcs
  attr(b_est, "nfpc") <- nfpc

  dir.create(file.path(results_wd, setting, scen, "b_tru"),
             showWarnings = FALSE)
  saveRDS(b_est, file = file.path(results_wd, setting, scen, "b_tru",
                                  paste0("b", i, ".rds")))
  NULL

}

sim_scen(n = 100, k = 5, r = 0.9, nk = 0.9, ce = 2, i = 1)
sim_scen(n = 100, k = 15, r = 0.9, nk = 0.9, ce = 2, i = 1)
