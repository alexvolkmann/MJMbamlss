#' Simulation Helper Function - Evaluate the Simulation for JMbamlss Setting
#'
#' This function evaluates the results for a given folder of JMbamlss model
#' fits.
#'
#' @param wd Path to simulations folder.
#' @param model_wd Simulation setting folder where the models are saved.
#' @param data_wd Simulation data folder.
#' @param name Name for description of the simulation setting.
#' @param rds Objects are saved as .rds files (for backwards compatibility when
#'   .Rdata files were used). Defaults to TRUE.
sim_jmbamlss_eval <- function(wd, model_wd, data_wd, name, rds = TRUE) {

  models <- list.files(path = paste0(wd, model_wd))
  if (rds) {
    models <- models[grep("\\.rds", models)]
  } else {
    rmv <- grep("\\.rds", models)
    if(length(rmv) > 0) {
      models <- models[-grep("\\.rds", models)]
    }
  }
  list_to_compare <- sim_bamlss_predict(models, wd, model_wd, data_wd, rds)

  it_list <- sim_results(lapply(list_to_compare, "[[", "predictions"),
                         lapply(list_to_compare, "[[", "simulations"),
                         name = name)
  do.call(rbind, Map(cbind, it = sub("\\.Rdata", "", names(it_list)), it_list))

}

#' Simulation Helper Function - Evaluate the Simulation for JMbayes Setting
#'
#' This function evaluates the results for a given folder of JMbayes model
#' fits.
#'
#' @param wd Path to simulations folder.
#' @param model_wd Simulation setting folder where the models are saved.
#' @param data_wd Simulation data folder.
#' @param name Name for description of the simulation setting.
#' @param rds Objects are saved as .rds files (for backwards compatibility when
#'   .Rdata files were used). Defaults to TRUE.
sim_jmbayes_eval <- function(wd, model_wd, data_wd, name, rds = TRUE) {

  models <- list.files(path = paste0(wd, model_wd))
  if (rds) {
    models <- models[grep("\\.rds", models)]
  } else {
    models <- models[-grep("\\.rds", models)]
  }
  list_to_compare <- sim_jmb_predict(models, wd, model_wd, data_wd, rds)

  it_list <- sim_results(lapply(list_to_compare, "[[", "predictions"),
                         lapply(list_to_compare, "[[", "simulations"),
                         name = name)
  do.call(rbind, Map(cbind, it = sub("\\.Rdata", "", names(it_list)), it_list))

}

sim_results <- function(result_list, dat_list, name) {

  n_dim <- length(levels(dat_list[[1]]$mu$marker))
  mapply(function (est, sim) {

    # Bias, rMSE, Coverage
    eval_lambga <- data.frame(
      type = c("Bias", "rMSE", "Coverage"),
      model = name,
      predictor = "lambga",
      marker = "all",
      t = "all",
      value = c(mean(est$lambga$Mean - sim$lambga),
                sqrt(mean((est$lambga$Mean - sim$lambga)^2)),
                mean(est$lambga[, 1] < sim$lambga &
                       est$lambga[, 3] > sim$lambga)))

    eval_alpha = data.frame(
      type = rep(c("Bias", "rMSE", "Coverage"), each = n_dim),
      model = name,
      predictor = "alpha",
      marker = paste0("m", seq_len(n_dim)),
      t = "all",
      value = c(est$alpha$Mean - sim$alpha,
                sqrt((est$alpha$Mean - sim$alpha)^2),
                as.numeric(est$alpha[, 1] < sim$alpha &
                             est$alpha[, 3] > sim$alpha)))

    eval_sigma = data.frame(
      type = rep(c("Bias", "rMSE", "Coverage"), each = n_dim),
      model = name,
      predictor = "sigma",
      marker = paste0("m", seq_len(n_dim)),
      t = "all",
      value = c(est$sigma$Mean - sim$sigma,
                sqrt((est$sigma$Mean - sim$sigma)^2),
                as.numeric(est$sigma[, 1] < sim$sigma &
                             est$sigma[, 3] > sim$sigma)))

    eval_mu <- data.frame(
      type = rep(c("Bias", "rMSE", "Coverage"), each = n_dim),
      model =  name,
      predictor = "mu",
      marker = paste0("m", seq_len(n_dim)),
      t = "all",
      value = c(mapply(function (e, s) {
                  mean(e - s)
                }, e = split(est$mu$Mean, est$mu$marker),
                s = split(sim$mu$mu, sim$mu$marker)),
                mapply(function (e, s) {
                  sqrt(mean((e - s)^2))
                }, e = split(est$mu$Mean, est$mu$marker),
                s = split(sim$mu$mu, sim$mu$marker)),
                mapply(function (l, u, s) {
                  mean(as.numeric(l < s & u > s))
                }, l = split(est$mu[, 1], est$mu$marker),
                u = split(est$mu[, 3], est$mu$marker),
                s = split(sim$mu$mu, sim$mu$marker))))

    sim_marker <- split(sim$mu_long, est$mu_long$marker)
    eval_mu_long <- data.frame(
      type = rep(c("Bias", "rMSE", "Coverage"), each = n_dim*101),
      model =  name,
      predictor = "mu_long",
      marker = paste0("m", seq(n_dim)),
      t = rep(rep(seq(0, 1, by = 0.01), each = n_dim), times = 3),
      value = c(c(sapply(seq(0, 1, by = 0.01), function (t) {
                  mapply(function(e, s) {
                    same_t <- e$obstime == t
                    mean(e$Mean[same_t] - s[same_t])
                  }, e = split(est$mu_long, est$mu_long$marker), s = sim_marker)
                })),
                c(sapply(seq(0, 1, by = 0.01), function (t) {
                  mapply(function(e, s) {
                    same_t <- e$obstime == t
                    sqrt(mean((e$Mean[same_t] - s[same_t])^2))
                  }, e = split(est$mu_long, est$mu_long$marker), s = sim_marker)
                })),
                c(sapply(seq(0, 1, by = 0.01), function (t) {
                  mapply(function(e, s) {
                    same_t <- e$obstime == t
                    mean(e[same_t, 1] < s[same_t] & e[same_t, 3] > s[same_t])
                  }, e = split(est$mu_long, est$mu_long$marker), s = sim_marker)
                }))))

    eval_lambga_long <- data.frame(
      type =  rep(c("Bias", "rMSE", "Coverage"), each = 101),
      model = name,
      predictor = "lambga_long",
      marker = "all",
      t = rep(seq(0, 1, by = 0.01), times = 3),
      value = c(c(sapply(seq(0, 1, by = 0.01), function (t) {
        e <- split(est$lambga_long, est$mu_long$marker)[[1]]
        # Following line should be replaced
        e$obstime <- split(est$mu_long, est$mu_long$marker)[[1]]$obstime
        s <- sim_marker[[1]]
        same_t <- e$obstime == t
        mean(e$Mean[same_t] - s[same_t])
      })),
      c(sapply(seq(0, 1, by = 0.01), function (t) {
        e <- split(est$lambga_long, est$mu_long$marker)[[1]]
        # Following line should be replaced
        e$obstime <- split(est$mu_long, est$mu_long$marker)[[1]]$obstime
        s <- sim_marker[[1]]
        same_t <- e$obstime == t
        sqrt(mean((e$Mean[same_t] - s[same_t])^2))
      })),
      c(sapply(seq(0, 1, by = 0.01), function (t) {
        e <- split(est$lambga_long, est$mu_long$marker)[[1]]
        # Following line should be replaced
        e$obstime <- split(est$mu_long, est$mu_long$marker)[[1]]$obstime
        s <- sim_marker[[1]]
        same_t <- e$obstime == t
        mean(e[same_t, 1] < s[same_t] & e[same_t, 3] > s[same_t])
      }))))

    eval_lambga_event <- data.frame(
      type = c("Bias", "rMSE", "Coverage"),
      model = name,
      predictor = "lambga_event",
      marker = "all",
      t = "all",
      value = c(mean(est$lambga_event$Mean - sim$lambga_event),
                sqrt(mean((est$lambga_event$Mean - sim$lambga_event)^2)),
                mean(est$lambga_event[, 1] < sim$lambga_event &
                       est$lambga_event[, 3] > sim$lambga_event)))

    rbind(eval_lambga, eval_alpha, eval_sigma, eval_mu, eval_mu_long,
          eval_lambga_long, eval_lambga_event)

  }, est = result_list, sim = dat_list, SIMPLIFY = FALSE)

}


sim_bamlss_predict_i <- function(m, wd, model_wd, data_wd, rds = TRUE,
                                 old = FALSE) {


  # Load the data set and extract information about it
  if (rds) {
    b_est <- readRDS(paste0(wd, model_wd, m))
    d_rirs <- readRDS(paste0(wd, data_wd, "d", substr(m, 2, 4), ".rds"))
  } else {
    load(paste0(wd, model_wd, m))
    load(paste0(wd, data_wd, "d", substr(m, 2, 4), ".Rdata"))
  }
  nodupl_ids <- which(!duplicated(b_est$model.frame[, c("id", "obstime")]))
  marks <- which(!duplicated(b_est$model.frame$marker))
  ids <- which(!duplicated(b_est$model.frame[, "id"]))
  d_rirs_obs <- d_rirs$data
  d_rirs_obs$survtime <- d_rirs_obs$obstime
  d_rirs_long <- d_rirs$data_full
  d_rirs_long$survtime <- d_rirs_long$obstime

  # Handle standardized survival matrices
  if (old) {

    # Standardization constants
    long_bar <- unique(samples(b_est)[, grep("long_bar",
                                             colnames(samples(b_est)))])
    long_sds <- unique(samples(b_est)[, grep("long_sds",
                                             colnames(samples(b_est)))])

    # Standardization constants for gamma
    gamma_bar <- b_est$x$gamma$smooth.construct$model.matrix$w_bar
    gamma_sds <- b_est$x$gamma$smooth.construct$model.matrix$w_sd
    if (length(b_est$x$gamma$smooth.construct) > 1) {
      warning("Accounting for standardization of smooths not yet implemented.")
    }
    std_x <- attr(b_est$terms$gamma, "term.labels")
    if (length(std_x) > 1) {
      warning("Accounting for standardization of more than one covariate not",
              " yet implemented.")
    }

    # Multiply the samples with standardization factor
    col_gamma <- grep(paste0("gamma\\.p\\.model\\.matrix\\.", std_x),
                      colnames(samples(b_est)))
    b_est$samples[[1]][, col_gamma] <- b_est$samples[[1]][, col_gamma] /
      gamma_sds
    col_alpha <- grep("alpha\\.", colnames(samples(b_est)))[
      seq_len(length(marks))]
    for (i in seq_along(col_alpha)) {
      b_est$samples[[1]][, col_alpha[i]] <-
        b_est$samples[[1]][, col_alpha[i]] / long_sds[i]
    }

    # Calculate standardization adaptation for intercept in each sample
    mcmc_int <- samples(b_est)[, col_gamma] * gamma_bar +
      apply(samples(b_est)[, col_alpha], 1, "%*%", long_bar)

    col_gamma_int <- grep("gamma\\.p\\.model\\.matrix\\.\\(Intercept",
                          colnames(samples(b_est)))
    b_est$samples[[1]][, col_gamma_int] <-
      b_est$samples[[1]][, col_gamma_int] - mcmc_int

  }

  # Extract MCMC samples to calculate the survival part
  # All longitudinal observation points
  mcmc_lambda <- as.matrix(predict(b_est, model = "lambda",
                                   newdata = d_rirs_obs,
                                   FUN = function(x) {x})[nodupl_ids, ])
  mcmc_gamma <- as.matrix(predict(b_est, model="gamma",
                                  newdata = d_rirs_obs,
                                  FUN = function(x) {x})[nodupl_ids, ])
  mcmc_lambga <- mcmc_gamma + mcmc_lambda
  # Longitudinal grid
  # Should be improved here. Only the first marker in d_rirs_long is needed.
  mcmc_lambda_long <- as.matrix(predict(b_est, model = "lambda",
                                        newdata = d_rirs_long,
                                        FUN = function(x) {x}))
  mcmc_gamma_long <- as.matrix(predict(b_est, model="gamma",
                                       newdata = d_rirs_long,
                                       FUN = function(x) {x}))
  mcmc_lambga_long <- mcmc_gamma_long + mcmc_lambda_long
  # Individual event times
  mcmc_lambda_event <- as.matrix(predict(b_est, model = "lambda",
                                         FUN = function(x) {x})[ids, ])
  mcmc_gamma_event <- as.matrix(predict(b_est, model="gamma",
                                        FUN = function(x) {x})[ids, ])
  mcmc_lambga_event <- mcmc_gamma_event + mcmc_lambda_event

  # Output list
  list("predictions" =  list(
    "lambga" = data.frame("2.5%" = apply(mcmc_lambga, 1, quantile,
                                         probs = 0.025),
                          "Mean" = rowMeans(mcmc_lambga),
                          "97.5%" = apply(mcmc_lambga, 1, quantile,
                                          probs = 0.975)),
    "alpha" = cbind(predict(b_est, model = "alpha", FUN = c95)[marks, ],
                    data.frame("marker" = b_est$model.frame$marker[marks])),
    "mu" = cbind(predict(b_est, model = "mu", FUN = c95),
                 data.frame("marker" = b_est$model.frame$marker,
                            "obstime" = d_rirs$data$obstime,
                            "id" = d_rirs$data$id)),
    "sigma" = cbind(predict(b_est, model = "sigma", FUN = c95)[marks, ],
                    data.frame("marker" = b_est$model.frame$marker[marks])),
    "mu_long" = cbind(predict(b_est, model = "mu", FUN = c95,
                              newdata = d_rirs$data_full),
                      data.frame("marker" = d_rirs$data_full$marker,
                                 "obstime" = d_rirs$data_full$obstime,
                                 "id" = d_rirs$data_full$id)),
    "lambga_long" = data.frame("2.5%" = apply(mcmc_lambga_long, 1, quantile,
                                              probs = 0.025),
                               "Mean" = rowMeans(mcmc_lambga_long),
                               "97.5%" = apply(mcmc_lambga_long, 1, quantile,
                                               probs = 0.975),
                               "marker" = d_rirs$data_full$marker,
                               "obstime" = d_rirs$data_full$obstime,
                               "id" = d_rirs$data_full$id),
    "lambga_event" = data.frame("2.5%" = apply(mcmc_lambga_event, 1, quantile,
                                         probs = 0.025),
                          "Mean" = rowMeans(mcmc_lambga_event),
                          "97.5%" = apply(mcmc_lambga_event, 1, quantile,
                                          probs = 0.975))
  ),
  "simulations" = list(
    "lambga" = rowSums(d_rirs$data[nodupl_ids, c("lambda", "gamma")]),
    "alpha" = d_rirs$data$alpha[marks],
    "mu" = d_rirs$data[, c("mu", "marker")],
    "sigma" = d_rirs$data$sigma[marks],
    "mu_long" = d_rirs$data_full$mu,
    "lambga_long" = data.frame(
      "lambga" = rowSums(d_rirs_long[, c("lambda", "gamma")]),
      "t" = d_rirs_long$survtime),
    "lambga_event" = rowSums(d_rirs$data[ids, c("lambda", "gamma")])
  ))

}

#' Simulation Helper Function - Predict the Results for bamlss-Models
#'
#' This function takes all the models listed in a folder and predicts the fit.
#'
#' @param m Vector containing all the file names of the models.
#' @param wd Path to simulations folder.
#' @param model_wd Simulation setting folder where the models are saved.
#' @param data_wd Simulation data folder.
#' @param rds Objects are saved as .rds files (for backwards compatibility when
#'   .Rdata files were used). Defaults to TRUE.
#' @param old Simulated data sets before Version 0.0.3 (samples need to be adapted
#'   for standardized survival matrices). Defaults to FALSE.
sim_bamlss_predict <- Vectorize(sim_bamlss_predict_i, vectorize.args = "m",
                                SIMPLIFY = FALSE)


sim_jmb_predict_i <- function(m, wd, model_wd, data_wd, rds = TRUE,
                              gamma_timeconst = TRUE) {

  # To not get a CMD check note for undefined global function
  # Remove if problematic
  survtime <- NULL

  # Load the fitted model and the original data
  if (rds) {
    jmb <- readRDS(paste0(wd, model_wd, m))
    d_rirs <- readRDS(paste0(wd, data_wd, "d", substr(m, 4, 6), ".rds"))
  } else {
    load(paste0(wd, model_wd, m))
    load(paste0(wd, data_wd, "d", substr(m, 4, 6), ".Rdata"))
  }

  # Problem if exogenous covariates in gamma: X_long for survival is not correct
  if(!gamma_timeconst) warning("Time constant gamma assumed")

  nodupl_ids <- which(!duplicated(d_rirs$data[, c("id", "obstime")]))
  n_dim <- length(levels(d_rirs$data$marker))
  marks <- which(!duplicated(d_rirs$data$marker))
  ids <- which(!duplicated(d_rirs$data[, "id"]))

  # Longitudinal fits
  X <- jmb$model_data$X
  Z <- jmb$model_data$Z
  B <- jmb$mcmc$b[[1]]
  n_re <- ncol(Z[[1]])
  mcmc_mu <- do.call(rbind, lapply(seq_len(n_dim), function (dim) {
    tcrossprod(X[[dim]], jmb$mcmc[[paste0("betas", dim)]][[1]]) +
      t(sapply(seq_len(nrow(Z[[dim]])), function (i) {
        Z[[dim]][i, ] %*% B[jmb$model_data$idL[[dim]][i],
                            (dim - 1)*n_re + seq_len(n_re), ]
      }))
  }))

  X_long <- split.data.frame(
    stats::model.matrix(formula(jmb$model_info$terms$terms_FE[[1]])[-2],
                        data = d_rirs$data_full),
    d_rirs$data_full$marker)
  Z_long <- split.data.frame(
    stats::model.matrix(formula(jmb$model_info$terms$terms_RE[[1]]),
                        data = d_rirs$data_full),
    d_rirs$data_full$marker)
  id_long <- split(d_rirs$data_full$id, d_rirs$data_full$marker)
  mcmc_mu_long <- do.call(rbind, lapply(seq_len(n_dim), function (dim) {
    tcrossprod(X_long[[dim]], jmb$mcmc[[paste0("betas", dim)]][[1]]) +
      t(sapply(seq_len(nrow(Z_long[[dim]])), function (i) {
        Z_long[[dim]][i, ] %*% B[id_long[[dim]][i],
                                 (dim - 1)*n_re + seq_len(n_re), ]
      }))
  }))

  # Survival fits
  kn <- mgcv::smoothCon(mgcv::s(survtime, k = 20, bs = "ps"),
                        data = d_rirs$data)[[1]]$knots
  Z <- splines::splineDesign(knots = kn, x = d_rirs$data$obstime, ord = 4,
                             outer.ok = TRUE)
  Z_long <- splines::splineDesign(knots = kn, x = d_rirs$data_full$obstime,
                                  ord = 4, outer.ok = TRUE)
  Z_event <- splines::splineDesign(knots = kn, x = d_rirs$data$survtime[ids],
                                   ord = 4, outer.ok = TRUE)
  X <- jmb$model_data$W_h[unlist(jmb$model_data$idL), , drop = FALSE]
  X_long <- jmb$model_data$W_h[as.numeric(d_rirs$data_full$id), , drop = FALSE]
  X_event <- jmb$model_data$W_h
  B <- jmb$mcmc$bs_gammas[[1]]
  Beta <- jmb$mcmc$gammas[[1]]
  mcmc_lambga <- (tcrossprod(Z, B) + tcrossprod(X, Beta))[nodupl_ids, ]
  mcmc_lambga_long <- (tcrossprod(Z_long, B) + tcrossprod(X_long, Beta))
  mcmc_lambga_event <- (tcrossprod(Z_event, B) + tcrossprod(X_event, Beta))

  list("predictions" = list(
        "lambga" = data.frame("2.5%" = apply(mcmc_lambga, 1, quantile,
                                             probs = 0.025),
                              "Mean" = rowMeans(mcmc_lambga),
                              "97.5%" = apply(mcmc_lambga, 1, quantile,
                                              probs = 0.975)),
        "alpha" = data.frame("2.5%" = jmb$statistics$CI_low$alphas,
                             "Mean" = jmb$statistics$Mean$alphas,
                             "97.5%" = jmb$statistics$CI_upp$alphas,
                             "marker" = factor(paste0("m", seq_len(n_dim)))),
        "mu" = data.frame("2.5%" = apply(mcmc_mu, 1, quantile,
                                         probs = 0.025),
                          "Mean" = rowMeans(mcmc_mu),
                          "97.5%" = apply(mcmc_mu, 1, quantile,
                                          probs = 0.975),
                          "marker" = d_rirs$data$marker,
                          "obstime" = d_rirs$data$obstime,
                          "id" = d_rirs$data$id),
        "sigma" = data.frame("2.5%" = log(jmb$statistics$CI_low$sigmas),
                             "Mean" = log(jmb$statistics$Mean$sigmas),
                             "97.5%" = log(jmb$statistics$CI_upp$sigmas),
                             "marker" = factor(paste0("m", seq_len(n_dim)))),
        "mu_long" = data.frame("2.5%" = apply(mcmc_mu_long, 1, quantile,
                                              probs = 0.025),
                               "Mean" = rowMeans(mcmc_mu_long),
                               "97.5%" = apply(mcmc_mu_long, 1, quantile,
                                               probs = 0.975),
                               "marker" = d_rirs$data_full$marker,
                               "obstime" = d_rirs$data_full$obstime,
                               "id" = d_rirs$data_full$id),
        "lambga_long" = data.frame("2.5%" = apply(mcmc_lambga_long, 1, quantile,
                                                  probs = 0.025),
                                   "Mean" = rowMeans(mcmc_lambga_long),
                                   "97.5%" = apply(mcmc_lambga_long, 1,
                                                   quantile, probs = 0.975)),
        "lambga_event" = data.frame("2.5%" = apply(mcmc_lambga_event, 1,
                                                   quantile, probs = 0.025),
                                    "Mean" = rowMeans(mcmc_lambga_event),
                                    "97.5%" = apply(mcmc_lambga_event, 1,
                                                    quantile, probs = 0.975))),
       "simulations" = list(
         "lambga" = rowSums(d_rirs$data[nodupl_ids, c("lambda", "gamma")]),
         "alpha" = d_rirs$data$alpha[marks],
         "mu" = d_rirs$data[, c("mu", "marker")],
         "sigma" = d_rirs$data$sigma[marks],
         "mu_long" = d_rirs$data_full$mu,
         "lambga_long" = data.frame(
           "lambga" = rowSums(d_rirs$data_full[, c("lambda", "gamma")]),
           "t" = d_rirs$data_full$obstime),
         "lambga_event" = rowSums(d_rirs$data[ids, c("lambda", "gamma")])
       ))


}

#' Simulation Helper Function - Predict the Results for JMbayes-Models
#'
#' This function takes all the models listed in a folder and predicts the fit.
#'
#' @param m Vector containing all the file names of the models.
#' @param wd Path to simulations folder.
#' @param model_wd Simulation setting folder where the models are saved.
#' @param data_wd Simulation data folder.
#' @param rds Objects are saved as .rds files (for backwards compatibility when
#'   .Rdata files were used). Defaults to TRUE.
#' @param gamma_timeconst Only implemented for timeconstant gamma predictors. If
#'   FALSE a warning message is returned.
sim_jmb_predict <- Vectorize(sim_jmb_predict_i, vectorize.args = "m",
                                SIMPLIFY = FALSE)
