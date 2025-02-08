#' Function to simulate data and estimate the model fits
#'
#' The function is hard coded so that only multiples of 3 can be used in as K.
#' This is connected to the specification of the types argument, which uses 3
#' "Own2" basis functions.
#'
#' @param it Number of iteration in the simulation scenario.
#' @param n Number of subjects in the simulated data set.
#' @param k Number of longitudinal markers. At the current implementation, can
#'  only be a multiple of 3.
#' @param r Correlation used to construct the Toeplitz matrix used for the
#'  covariance of the univariate basis functions.
#' @param nk Probability of missingness of the longitudinal observations.
#'  Corresponds to the argument probmiss of function MJMbamlss::simMultiJM().
#' @param ce Factor changing the uniform censoring interval. Corresponds to the
#'  argument maxfac of function MJMbamlss::simMultiJM().
#' @param type Type of univariate functional principal components to simulate.
#'  Available options are the types available to the funData::eFun() function
#'  (e.g. "Fourier", "Poly", or "Wiener"). Additionally, the two sets of basis
#'  functions used in simulation scenario II ("Own1" and "Own2" for the
#'  different dimensions) are possible. Defaults to "Own2".
#' @param nb Number of univariate basis functions per dimension in the
#'  construction of the multivariate FPCs. Defaults to 3.
#' @param alpha Vector of length 3 specifying the values of the linear predictor
#'  of alpha.
#' @param gamma Intercept for the baseline hazard. Defaults to a function
#'  depending on the number of longitudinal markers k.
#' @param mfpc_true TRUE if the true, data generating MFPC basis should be used.
#'  If FALSE, the MFPC basis is estimated. Defaults to TRUE.
#' @param mfpc_truncation Truncation order for the MFPC basis. Numeric values
#'  like 0.95 truncate after 95% explained variation and 1L includes the entire
#'  MFPC basis.
#' @param JMBayes TRUE if the package JMBayes2 is to be used for the model fit.
#'  Defaults to FALSE, which corresponds to MJMbamlss using MFPCs.
#' @param results_wd String specifying where the folders for saving the objects
#'   should be created.
sim_fun <- function (it, n, k, r, nk, ce, type = "Own2", nb = 3, alpha,
                     gamma = -2 - (k %/% 3), mfpc_true = TRUE,
                     mfpc_truncation = 1L, JMBayes = FALSE, results_wd) {

  require(MJMbamlss)
  require(funData)
  require(tidyverse)
  if (any(k %% 3 != 0, length(alpha) %% 3 != 0)) {
    stop("Arguments k needs to be a multiple of 3 and alpha needs to ",
         "be of length 3!")
  }

  # Name for scenario and model
  scen <- paste("set", type, n, k, r, nk, ce, nb, paste0(alpha, collapse = "_"),
                gamma, sep = "_")
  model <- paste("mod", if(JMBayes) "jmb" else "mjm", mfpc_true,
                 mfpc_truncation,sep = "_")

  # Create the folder and subfolders
  dir.create(file.path(results_wd, scen),
             showWarnings = FALSE)
  dir.create(file.path(results_wd, scen, "data"),
             showWarnings = FALSE)
  dir.create(file.path(results_wd, scen, model),
             showWarnings = FALSE)

  # Generate (if it doesn't yet exist) or load the data
  set.seed(it)
  data_file <- file.path(results_wd, scen, "data", paste0("d", it, ".rds"))
  if (!file.exists(data_file)) {
    d_sim <- generate_data(results_wd, scen, k, type, nb, r, it, n, nk, ce,
                            gamma, alpha)
  } else {
    d_sim <- readRDS(data_file)
  }
  d <- d_sim$data


  # MODEL ESTIMATION
  if (JMBayes) {

    # JMBayes2 MODEL ESTIMATION
    require(JMbayes2)

    # Get the quantile-based knots for comparability
    kn <- mgcv::smoothCon(mgcv::s(survtime, k = 20, bs = "ps"),
                          data = d)[[1]]$knots

    # Estimate the model using JMbayes
    d_sim_jmb <- d %>%
      pivot_wider(names_from = marker, values_from = y)
    t_jm <- system.time({
      CoxFit <- coxph(Surv(survtime, event) ~ x3,
                      data = d_sim$data_short %>% filter(marker == "m1"))
      lme_list <- lapply(seq_len(k), function (i) {
          form <- as.formula(paste0(paste0("m", i), "~ obstime * x3"))
          assign(
            paste0("lm", i),
            lme(fixed = form,
                random = ~ bs(obstime, df = 3, Boundary.knots = c(0, 1)) | id,
                data = d_sim_jmb, na.action = na.omit,
                control = lmeControl(opt = "optim")))
      })
      jmb <- jm(CoxFit, lme_list, time_var = "obstime",
                n_iter = 5500L, n_burnin = 500L, n_thin = 5L,
                cores = 1, n_chains = 1,
                GK_k = 7, Bsplines_degree = 3, diff = 2, knots = list(kn),
                save_random_effects = TRUE)
    })
    attr(jmb, "comp_time") <- t_jm
    saveRDS(jmb, file = file.path(results_wd, scen, model,
                                  paste0("jmb", it, ".rds")))

  } else {

    # MJM MODEL ESTIMATION
    t_est <- system.time({
      # Use the true data generating FPCs or estimate them
      if (mfpc_true) {
        mfpcs <- readRDS(file.path(results_wd, scen, "mfpc.rds"))
      } else {
        mfpcs <- estimate_mfpc_basis(data = d, n_min = 3, minmaxobs_long = 0.1,
                                     npcs = 3, n_basis = 10)
      }

      # Truncate the basis if needed
      if (mfpc_truncation == 1L) {
        nfpc <- length(which(mfpcs$values > 0))
        mfpc_list <- lapply(seq_len(nfpc), function (it) {
          list(functions = extractObs(mfpcs$functions, it),
               values = mfpcs$values[it])
        })
      } else {
        vals <- which(mfpcs$values > 0)
        nfpc <- min(which(
          cumsum(mfpcs$values[vals]) / sum(mfpcs$values[vals]) >
            mfpc_truncation))
        mfpc_list <- lapply(vals[seq_len(nfpc)], function (i, mfpca = mfpcs){
          list(functions = extractObs(mfpcs$functions, i),
               values = mfpcs$values[i])
        })
      }

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
                   "mfpc_list[[", x, "]], 'scale' = FALSE))")
          }), collapse = " + "))),
        sigma ~ -1 + marker,
        alpha ~ -1 + marker
      )

      # Run the model
      set.seed(it)

      b_est <- bamlss(form, family = mjm_bamlss,
                      data = d,
                      timevar = "obstime", maxit = 1500, n.iter = 5500,
                      burnin = 500, thin = 5)
    })
    attr(b_est, "comp_time") <- t_est
    attr(b_est, "FPCs") <- mfpcs
    attr(b_est, "nfpc") <- nfpc

    saveRDS(b_est, file = file.path(results_wd, scen, model,
                                    paste0("b", it, ".rds")))
  }

  NULL

}


# Data Simulation ---------------------------------------------------------

# Function to generate the data
generate_data <- function (results_wd, scen, k, type, nb, r, it, n, nk, ce,
                           gamma, alpha) {

  # Functional observation grid
  argvals <- seq(0, 1, by = 0.01)

  # Create the MFPCs if it does not exist yet
  mfpc_file <- file.path(results_wd, scen, "mfpc.rds")
  if (!file.exists(mfpc_file)) {
    mfpca <- create_mfpca(argvals, k, type, nb, r, mfpc_file)
  } else {
    mfpca <- readRDS(mfpc_file)
  }

  # Simulate the data
  set.seed(it)
  d_sim <- simMultiJM(nsub = n, times = argvals,
                      max_obs = 15, probmiss = nk, maxfac = ce,
                      nmark = k, long_assoc = "FPC", M = k * nb,
                      FPC_bases = mfpca$functions, FPC_evals = mfpca$values,
                      ncovar = 2,
                      lambda = function(t, x) {
                        1.65 * t^(0.65)
                      },
                      gamma = function(x) {
                        gamma + 0.3*x[, 3]
                      },
                      alpha = rep(list(
                        function(t, x) alpha[1],
                        function(t, x) alpha[2],
                        function(t, x) alpha[3]), k %/% 3),
                      mu = rep(list(function(t, x, r){
                        0 + 1*t + 0.3*x[, 3] + 0.3*t*x[, 3]
                      }), k),
                      sigma = function(t, x) {
                        log(0.06) + 0*t
                      },
                      tmax = NULL, seed = NULL,
                      full = TRUE, file = NULL)

  saveRDS(d_sim,
          file = file.path(results_wd, scen, "data",
                           paste0("d", it, ".rds")))
  d_sim

}

# Function to generate the MFPCA object for data generation
create_mfpca <- function (argvals, k, type, nb, r, mfpc_file) {

  # Type of basis functions to use for data simulation
  if (type == "Own1") {

    # Same eigenfunctions as in simulation scenario II, marker 1
    x <- seq(0, 1, by = 0.1)
    m1sp1 <- splinefun(
      x, c(0.3, 0.5, 0.6, 0.5, 0.4, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5))
    m1sp2 <- splinefun(
      x, c(0, 0, 0, 0.3, 0.5, 0.6, 0.5, 0.4, 0.5, 0.5, 0.5))
    m1sp3 <- splinefun(
      x, c(0, 0, 0, 0, 0, 0, 0.3, 0.5, 0.6, 0.5, 0.4))
    m1 <- funData(
      argvals = argvals,
      X = matrix(c(m1sp1(argvals), m1sp2(argvals), m1sp3(argvals)),
                 nrow = 3, byrow = TRUE))

    B <- rep(list(m1), k)

  } else if (type == "Own2") {

    # Simulation scenario II, marker2
    x <- seq(0, 1, by = 0.1)
    m2sp1 <- splinefun(x,
                       c(0.5, 1, 0.7, 0.5, 0.3, 0.2, 0.1, 0.05, 0, 0, 0))
    m2sp2 <- splinefun(x,
                       c(0, 0, 0, 0.5, 1, 0.7, 0.5, 0.3, 0.2, 0.1, 0.05))
    m2sp3 <- splinefun(x,
                       c(0, 0, 0, 0, 0, 0, 0.5, 1, 0.7, 0.5, 0.3))

    m2 <- funData(
      argvals = argvals,
      X = matrix(c(m2sp1(argvals), m2sp2(argvals), m2sp3(argvals)),
                 nrow = 3, byrow = TRUE))

    B <- rep(list(m2), k)

  } else {
    B <- rep(list(funData::eFun(argvals, nb, type = type)), k)
  }

  # Toeplitz covariance for basis coefficients
  COV <- toeplitz(toeplitzfun(1:(k*nb), rho = r))

  # Calculate MFPCs
  mfpca <- MFPCA_cov(cov = COV, basis_funs = B)
  mfpca$COV <- COV

  saveRDS(mfpca, file = mfpc_file)
  mfpca

}



# Prepare MFPC Basis for Model Fit ----------------------------------------

#' Function to estimate the MFPC basis
#'
#' @param data Simulated data to estimate the MFPC basis on.
#' @param n_min Minimum number of longitudinal observations for the MFPC
#'  estimation. Defaults to 3.
#' @param minmaxobs_long Minimum for the last (maximal) longitudinal time
#'  observation to be considered for the MFPC estimation. Defaults to 0.1.
#' @param npcs Number of univariate principal components to use in estimation.
#'  See preproc_MFPCA function, argument npc. Defaults to 3.
#' @param n_basis Number of b-splines basis functions in estimation. See
#'  preproc_MFPCA function, arguemnt nbasis. Defaults to 10.
#' @returns Estimated MFPCAfit object as outcome from function preproc_MFPCA.
estimate_mfpc_basis <- function (data, n_min = 3, minmaxobs_long = 0.1,
                                 npcs = 3, n_basis = 10) {

    # Remove observation with less than 3 longitudinal observations
    few_obs <- apply(table(data$id, data$marker), 1,
                     function (x) any(x < n_min))
    long_obs <- data %>%
      group_by(id, marker) %>%
      summarize(maxobs = max(obstime), .groups = "drop_last") %>%
      ungroup(marker) %>%
      summarize(minmaxobs = min(maxobs), .groups = "drop_last") %>%
      ungroup() %>%
      filter(minmaxobs > minmaxobs_long) %>%
      select(id) %>%
      unlist() %>%
      paste()
    take <- intersect(long_obs, paste(which(!few_obs)))

    mfpca_est <- preproc_MFPCA(
      data %>%
        filter(id %in% take) %>%
        droplevels(),
      uni_mean = "y ~ 1 + obstime + x3 + obstime:x3",
      npc = npcs, nbasis = n_basis
    )
    mfpca_est

}


# Function for Covariance Matrix ------------------------------------------

# Function to calculate the elements of the Toeplitz-Matrix
toeplitzfun <- Vectorize(function (x, rho) {
  exp(1 / (rho)^2 * (1/(20*3) - x/(20*3)))
}, vectorize.args = "x")



# Functions for Evaluation ------------------------------------------------

#' Function to Simplify the Evaluation of Simulation Scenarios
#'
#' This function is a helperfunction to evaluate all the predicted models of one
#' specific model specification and saves it in the respective folder.
#'
#' @param wd Path to simulations folder.
#' @param model_wd Simulation setting folder where the models are saved.
#' @param JMBayes TRUE if the package JMBayes2 is to be used for the model fit.
#'  Defaults to FALSE, which corresponds to MJMbamlss using MFPCs.
#'
#' @returns NULL. Saves the predictions and evaluations of the simulation runs.
evaluate_sim_setting <- function(wd, model_wd, JMBayes = FALSE) {


  m_wd <- if (substr(model_wd, nchar(model_wd)-1, nchar(model_wd)) != "\\") {
    file.path(model_wd, "")
  } else model_wd
  s_wd <- if (substr(wd, nchar(wd)-1, nchar(wd)) != "\\") {
    file.path(wd, "")
  } else wd

  # Create the folder for saving the evaluations
  dir.create(file.path(wd, "eval"), showWarnings = FALSE)
  dir.create(file.path(wd, "eval", "predictions"), showWarnings = FALSE)
  dir.create(file.path(wd, "eval", "evaluations"), showWarnings = FALSE)

  # Predict all of the models in the folder
  model_files <- list.files(path = file.path(wd, model_wd))

  if (JMBayes) {
    preds <- MJMbamlss:::sim_jmb_predict(
      m = model_files, wd = s_wd, model_wd = m_wd,
      data_wd = file.path("data", "")
    )
  } else {
    preds <- MJMbamlss:::sim_bamlss_predict(
      m = model_files, wd = s_wd, model_wd = m_wd,
      data_wd = file.path("data", "")
    )
  }
  saveRDS(preds, file = file.path(wd, "eval", "predictions",
                                  paste0("p_", model_wd,".rds")))

  # Evaluate all of the models
  it_list <- MJMbamlss:::sim_results(lapply(preds, "[[", "predictions"),
                                     lapply(preds, "[[", "simulations"),
                                     name = model_wd)
  evals <- do.call(rbind, Map(cbind, it = sub("\\.rds", "", names(it_list)),
                              it_list))
  saveRDS(evals,
          file = file.path(wd, "eval", "evaluations",
                           paste0("e_", model_wd, ".rds")))
  NULL

}

