
# Set Up Simulation: Data -------------------------------------------------



# Set up R session --------------------------------------------------------



# Specify location
# setwd()
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
# Uncomment to create the structure of the simulated objects
# dir.create(file.path(results_wd, setting), showWarnings = FALSE)
# dir.create(file.path(results_wd, setting, "data"), showWarnings = FALSE)
# dir.create(file.path(results_wd, setting,
#                      "/bamlss_tru"), showWarnings = FALSE)
# dir.create(file.path(results_wd, setting,
#                      "/bamlss_est1"), showWarnings = FALSE)
# dir.create(file.path(results_wd, setting, "/bamlss_est99"),
#            showWarnings = FALSE)
# dir.create(file.path(results_wd, setting, "/jmb"), showWarnings = FALSE)
# dir.create(file.path(results_wd, setting, "data", "mfpcs"),
#            showWarnings = FALSE)



# Different Specifications ------------------------------------------------

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
                results_wd, setting, "data", "mfpcs", paste0(
                  paste("mfpc", k, r, sep = "_"), ".rds"
                )))
  }
}


# Create the directories to save the data ---------------------------------

for (n in N) {
  for (k in K) {
    for (r in R) {
      for (nk in NK) {
        for (ce in CE) {
          dir.create(file.path(results_wd, setting, "data",
                               paste("d", n, k, r, nk, ce, sep = "_")),
                     showWarnings = FALSE)
        }
      }
    }
  }
}

# Simulation function -----------------------------------------------------
parallel_data <- function (i) {

  for (n in N) {
    for (k in K) {
      for (r in R) {
        for (nk in NK) {
          for (ce in CE) {

            # Load the corresponding MFPCs
            m <- readRDS(file.path(
              results_wd, setting, "data", "mfpcs", paste0(
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
                    file = file.path(results_wd, setting, "data",
                                     paste("d", n, k, r, nk, ce, sep = "_"),
                                     paste0("d", i, ".rds")))
            NULL

          }
        }
      }
    }
  }

}

# d_sim <- readRDS(file.path(results_wd, setting, "data",
#                            paste("d", n, k, r, nk, ce, sep = "_"),
#                            paste0("d", i, ".rds")))
# plot(survfit(Surv(survtime, event) ~ x3, data = d_sim$data_short[1:n, ]))
# table(d_sim$data_short$event[1:n])
# d_sim$data_short %>%
#   slice_head(n = n) %>%
#   mutate(cens = ifelse(survtime > 0.999 & event == 0, "full",
#                        ifelse(event == 0, "cens", "event"))) %>%
#   pull(cens) %>%
#   table()
# ggplot(data = d_sim$data_short %>% slice_head(n = n),
#        aes(x = id, y = survtime, colour = as.factor(event))) +
#   geom_point() +
#   guides(colour =guide_legend(title = "Event"))
# ggplot(d_sim$data %>% filter(marker %in% paste0("m", 1:5)),
#        aes(x = obstime, y = y, color = id)) +
#   geom_line() +
#   facet_grid(~marker) +
#   theme(legend.position = "none")
# d_sim$data %>%
#   group_by(marker, id) %>%
#   summarize(N = n()) %>%
#   pull(N) %>% median() #mean()
#   table()

# Simulation --------------------------------------------------------------

cat(paste("Parallelization:", start, "-", stop, "on", number_cores, "cores.\n"))
dat_sim <- mclapply(start:stop, parallel_data, mc.cores = number_cores)

