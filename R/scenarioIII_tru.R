
# Run Models with True MFPCs ----------------------------------------------


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

# Setting for the simulation
start <- 1
stop <- 10
number_cores <- 5
setting <- "scen_III"



# Different Specifications ------------------------------------------------

# Specifications are ordered so that estimation should get successively more
# computationally expensive as more data are available


# Number of individuals
N <- c(50, 250, 500)

# Numbers of longitudinal markers
K <- c(5, 10, 15)

# Correlation between markers
R <- c(0.9, 0.5, 0.2)

# Marker sparsity
NK <- c(0.9, 0.75, 0.5)

# Censoring rate
CE <- c(2, 3, 4)




# Function to estimate ----------------------------------------------------

bamlss_est <- function (n, k, r, nk, ce, i) {

  # Load the data
  d_sim <- readRDS(file.path(results_wd, setting, "data",
                             paste("d", n, k, r, nk, ce, sep = "_"),
                             paste0("d", i, ".rds")))
  d <- d_sim$data
  mfpcs <- d_sim$fpc_base

  # Use all true data generating FPCs
  nfpc <- 3#nObs(d_sim$fpc_base)
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
  saveRDS(b_est, file = file.path(results_wd, setting, "bamlss_tru",
                                  paste("d", n, k, r, nk, ce, sep = "_"),
                                  paste0("b", i, ".rds")))
  cat(i, "\\")

}



# Try out -----------------------------------------------------------------

bamlss_est(n = 50, k = 5, r = 0.9, nk = 0.9, ce = 2, i = 1)
