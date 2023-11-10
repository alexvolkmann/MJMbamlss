

# Simulation Scenario II - 22/12/09 ---------------------------------------


# JMbayes estimation


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



# Simulation function -----------------------------------------------------

parallel_bamlss_est <- function(i) {
  set.seed(i)
  
  # Load the data
  d_sim <- readRDS(file.path(results_wd, setting, "data",
                             paste0("d", i, ".rds")))
  
  # Get the quantile-based knots for comparability
  kn <- mgcv::smoothCon(mgcv::s(survtime, k = 20, bs = "ps"), 
                        data = d_sim$data)[[1]]$knots
  
  # Estimate the model using JMbayes
  d_sim_jmb <- d_sim$data %>% 
    pivot_wider(names_from = marker, values_from = y)
  try_obj <- try({
    t_jm <- system.time({
      CoxFit <- coxph(Surv(survtime, event) ~ x3, 
                      data = d_sim$data_short %>% filter(marker == "m1"))
      lm1 <- lme(m1 ~ obstime * x3, 
                 random = ~ bs(obstime, df = 3, Boundary.knots = c(0, 1)) | id, 
                 data = d_sim_jmb, na.action = na.omit,
                 control = lmeControl(opt = "optim"))
      lm2 <- lme(m2 ~ obstime * x3, 
                 random = ~ bs(obstime, df = 3, Boundary.knots = c(0, 1)) | id, 
                 data = d_sim_jmb, na.action = na.omit,
                 control = lmeControl(opt = "optim"))
      jmb <- jm(CoxFit, list(lm1, lm2), time_var = "obstime",
                n_iter = 5500L, n_burnin = 500L, n_thin = 5L, 
                cores = 1, n_chains = 1, 
                GK_k = 7, Bsplines_degree = 3, diff = 2, knots = list(kn),
                save_random_effects = TRUE)
    })
    attr(jmb, "comp_time") <- t_jm
    saveRDS(jmb, file = file.path(results_wd, setting, "jmb", 
                                    paste0("jmb", i, ".rds")))
      
  }, silent = TRUE)
  if(class(try_obj) == "try-error") try_obj else i
}



# Simulation --------------------------------------------------------------

cat(paste("Parallelization:", start, "-", stop, "on", number_cores, "cores.\n"))
simulation <- mclapply(start:stop, parallel_bamlss_est, mc.cores = number_cores)
