
# .libPaths(...)
# library(devtools)
# install_github("alexvolkmann/MJMbamlss")
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



# Simulation settings -----------------------------------------------------


# Specifications are ordered so that estimation should get successively more
# computationally expensive as more data are available


# Number of individuals
N <- c(150, 300, 500)

# Numbers of longitudinal markers
K <- c(3, 6, 9)

# Correlation between markers
R <- c(0.9, 0.5, 0.2)

# Marker sparsity
NK <- c(0.8, 0.6, 0.4)

# Censoring rate
CE <- c(1.5, 1.75, 2)

setting_args <- expand.grid(n = N, k = K, r = R, nk = NK, ce = CE)


# Outline
# Start with
# Full N, K
# R[c(1, 2)]
# NK[1], CE[1]


# Simulation Setting Function ---------------------------------------------

run_simulation_setting <- function(its, n, k, r, nk, ce, results_wd,
                                   gamma = -2 - (k %/% 3)) {

  # Run the models for all iterations
  for (i in its) {

    mjm_true <- MJMbamlss:::sim_fun(i, n, k, r, nk, ce, type = "Own2", nb = 3,
                        alpha = c(1, 0.75, 0.5),
                        gamma = -2 - (k %/% 3),
                        mfpc_true = TRUE,
                        mfpc_truncation = 1L,
                        JMBayes = FALSE,
                        results_wd = results_wd)

    jmb <- MJMbamlss:::sim_fun(i, n, k, r, nk, ce, type = "Own2", nb = 3,
                   alpha = c(1, 0.75, 0.5),
                   gamma = -2 - (k %/% 3),
                   mfpc_true = TRUE,
                   mfpc_truncation = 1L,
                   JMBayes = TRUE,
                   results_wd = results_wd)

    mjm_est <- MJMbamlss:::sim_fun(i, n, k, r, nk, ce, type = "Own2", nb = 3,
                       alpha = c(1, 0.75, 0.5),
                       gamma = -2 - (k %/% 3),
                       mfpc_true = FALSE,
                       mfpc_truncation = 1L,
                       JMBayes = FALSE,
                       results_wd = results_wd)

    mjm_tr99 <- MJMbamlss:::sim_fun(i, n, k, r, nk, ce, type = "Own2", nb = 3,
                        alpha = c(1, 0.75, 0.5),
                        mfpc_true = TRUE,
                        mfpc_truncation = 0.99,
                        JMBayes = FALSE,
                        results_wd = results_wd)

  }

  # Name for scenario and model
  scen <- paste("set", "Own2", n, k, r, nk, ce, nb,
                paste0(alpha, collapse = "_"),
                gamma, sep = "_")
  model <- paste("mod", if(JMBayes) "jmb" else "mjm", mfpc_true,
                 mfpc_truncation, sep = "_")
  eval_true <- MJMbamlss:::evaluate_sim_setting(
    wd = file.path(results_wd, scen),
    model_wd = paste("mod", "mjm", TRUE, 1, sep = "_"),
    JMBayes = FALSE
  )
  eval_jmb <- MJMbamlss:::evaluate_sim_setting(
    wd = file.path(results_wd, scen),
    model_wd = paste("mod", "jmb", TRUE, 1, sep = "_"),
    JMBayes = FALSE
  )
  eval_est <- MJMbamlss:::evaluate_sim_setting(
    wd = file.path(results_wd, scen),
    model_wd = paste("mod", "mjm", FALSE, 1, sep = "_"),
    JMBayes = FALSE
  )
  eval_est <- MJMbamlss:::evaluate_sim_setting(
    wd = file.path(results_wd, scen),
    model_wd = paste("mod", "mjm", FALSE, 0.99, sep = "_"),
    JMBayes = FALSE
  )

}



# Run Simulation Scenario -------------------------------------------------

# Arguments 'its' beginning with 100:109
run_simulation_setting(its = c(200, 201), n = 150, k = 3, r = 0.9,
                       nk = 0.3, ce = 2, results_wd = ".")
