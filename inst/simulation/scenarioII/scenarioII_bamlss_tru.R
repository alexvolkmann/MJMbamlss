

# Simulation Scenario II - 22/12/09 ---------------------------------------




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
number_cores <- 2
setting <- "scen_II"
Sys.time()
sessionInfo()



# Recreate true MFPCA -----------------------------------------------------

# Number of individuals and other quantities
argvals <- seq(0, 1, by = 0.01)
x <- seq(0, 1, by = 0.1)

# Random covariance matrix
# Set the eigenvalues but the eigenvectors are random
set.seed(1105)
p <- 6
P <- qr.Q(qr(matrix(rnorm(p^2), ncol = p)))
evals <- c(4, 3, 2, 1, 0.5, 0.2)
cov <- crossprod(P, P*(evals))


# Find spline functions
# Marker1
m1sp1 <- splinefun(x, c(0.3, 0.5, 0.6, 0.5, 0.4, 0.5, 0.5, 0.5, 0.5, 0.5, 
                        0.5))
m1sp2 <- splinefun(x, c(0, 0, 0, 0.3, 0.5, 0.6, 0.5, 0.4, 0.5, 0.5, 0.5))
m1sp3 <- splinefun(x, c(0, 0, 0, 0, 0, 0, 0.3, 0.5, 0.6, 0.5, 0.4))
# Marker2
m2sp1 <- splinefun(x, c(0.5, 1, 0.7, 0.5, 0.3, 0.2, 0.1, 0.05, 0, 0, 0))
m2sp2 <- splinefun(x, c(0, 0, 0, 0.5, 1, 0.7, 0.5, 0.3, 0.2, 0.1, 0.05))
m2sp3 <- splinefun(x, c(0, 0, 0, 0, 0, 0, 0.5, 1, 0.7, 0.5, 0.3))

m1 <- funData(argvals = argvals,
              X = matrix(c(m1sp1(argvals), m1sp2(argvals), m1sp3(argvals)),
                         nrow = 3, byrow = TRUE))
m2 <- funData(argvals = argvals,
              X = matrix(c(m2sp1(argvals), m2sp2(argvals), m2sp3(argvals)),
                         nrow = 3, byrow = TRUE))

# True multivariate covariance structure
m <- MFPCA_cov(cov = cov, basis_funs = list(m1, m2))



# Simulation function -----------------------------------------------------

parallel_bamlss_est <- function(i) {
  set.seed(i)
  
  # Load the data
  d_sim <- readRDS(file.path(results_wd, setting, "data",
                              paste0("d", i, ".rds")))
  
  try_obj <- try({
    
    # Use the true data generating FPCs
    # Use all FPCs 
    nfpc <- nObs(d_sim$fpc_base)
    mfpca_est_list <- lapply(seq_len(nfpc), function (i) {
      list(functions = extractObs(m$functions, i),
           values = m$values[i])
    })
    
    # Prepare objects for model fit
    f_est <- list(
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
    
    # Model fit
    t_est <- system.time(
      b_est <- bamlss(f_est, family = mjm_bamlss, 
                      data = d_sim$data, 
                      timevar = "obstime", maxit = 1500, n.iter = 5500,
                      burnin = 500, thin = 5)
    )
    attr(b_est, "comp_time") <- t_est
    attr(b_est, "FPCs") <- m
    attr(b_est, "nfpc") <- nfpc
    saveRDS(b_est, file = file.path(results_wd, setting, "bamlss_tru", 
                                    paste0("b", i, ".rds")))
  }, silent = TRUE)
  if(class(try_obj) == "try-error") try_obj else i
}



# Simulation --------------------------------------------------------------

cat(paste("Parallelization:", start, "-", stop, "on", number_cores, 
          "cores.\n"))
simulation <- mclapply(start:stop, parallel_bamlss_est, 
                       mc.cores = number_cores)
