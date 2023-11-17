
# Set Up Simulation: Data -------------------------------------------------



# Set up R session --------------------------------------------------------



# Specify location
# setwd()
# results_wd <-



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


# Setting for the simulation
start <- 100
stop <- 299
number_cores <- 5
setting <- "scen_I"
# Uncomment to create the structure of the simulated objects
# dir.create(file.path(results_wd, setting), showWarnings = FALSE)
# dir.create(file.path(results_wd, setting, "data"), showWarnings = FALSE)
# dir.create(file.path(results_wd, setting, "/bamlss_tru"), showWarnings = FALSE)
# dir.create(file.path(results_wd, setting, "/bamlss_est1"), showWarnings = FALSE)
# dir.create(file.path(results_wd, setting, "/bamlss_est99"),
#            showWarnings = FALSE)
# dir.create(file.path(results_wd, setting, "/jmb"), showWarnings = FALSE)
Sys.time()
sessionInfo()


# Objects for all clusters ------------------------------------------------

# Number of individuals
n <- 150

# Covariance matrix for the data generation
auto <- matrix(c(0.08, -0.07, -0.07, 0.9), ncol = 2)
cross <- matrix(rep(0.03, 4), ncol = 2)
cor <- matrix(c(0, 1, 0.75, 0.5, 0, 0,
                1, 0, 1, 0.75, 0.5, 0,
                0.75, 1, 0, 1, 0.75, 0.5,
                0.5, 0.75, 1, 0, 1, 0.75,
                0, 0.5, 0.75, 1, 0, 1,
                0, 0, 0.5, 0.75, 1, 0),
              ncol = 6)
cov <- kronecker(cor, cross) + kronecker(diag(c(1, 1.2, 1.4, 1.6, 1.8, 2)), auto)


# Simulation function -----------------------------------------------------

parallel_data <- function(i) {
  set.seed(i)

  # Simulate the data
  d_rirs <- simMultiJM(
    nsub = n, times = seq(0, 1, by = 0.01),
    max_obs = 15, probmiss = 0.75, maxfac = 1.75,
    nmark = 6, long_assoc = "param", M = NULL,
    FPC_bases = NULL, FPC_evals = NULL, mfpc_args = NULL,
    re_cov_mat = cov,
    ncovar = 2,
    lambda = function(t, x) {
      1.37 * t^(0.37)
    },
    gamma = function(x) {
      -1.5 + 0.48*x[, 3]
    },
    alpha = list(function(t, x) {
      1.5 + 0*t
    }, function(t, x) {
      0.6 + 0*t
    }, function(t, x) {
      0.3 + 0*t
    }, function(t, x) {
      -0.3 + 0*t
    }, function(t, x) {
      -0.6 + 0*t
    }, function(t, x) {
      -1.5 + 0*t
    }),
    mu = list(function(t, x, r){
      0 + 0.2*t - 0.25*x[, 3] - 0.05*t*x[, 3] +
        r[, 1] + r[, 2]*t
    }, function(t, x, r){
      0 + 0.2*t - 0.25*x[, 3] - 0.05*t*x[, 3] +
        r[, 3] + r[, 4]*t
    }, function(t, x, r){
      0 + 0.2*t - 0.25*x[, 3] - 0.05*t*x[, 3] +
        r[, 5] + r[, 6]*t
    }, function(t, x, r){
      0 + 0.2*t - 0.25*x[, 3] - 0.05*t*x[, 3] +
        r[, 7] + r[, 8]*t
    }, function(t, x, r){
      0 + 0.2*t - 0.25*x[, 3] - 0.05*t*x[, 3] +
        r[, 9] + r[, 10]*t
    }, function(t, x, r){
      0 + 0.2*t - 0.25*x[, 3] - 0.05*t*x[, 3] +
        r[, 11] + r[, 12]*t
    }),
    sigma = function(t, x) {
      log(0.06) + 0*t
    },
    tmax = NULL, seed = NULL,
    full = TRUE, file = NULL
  )

  # Uncomment to save the simulated data
  # saveRDS(d_rirs, file = file.path(results_wd, setting, "data",
  #                                  paste0("d", i, ".rds")))
  NULL
}



# Simulation --------------------------------------------------------------

cat(paste("Parallelization:", start, "-", stop, "on", number_cores, "cores.\n"))
mclapply(start:stop, parallel_data, mc.cores = number_cores)

