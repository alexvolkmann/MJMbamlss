

# Optimizer for MJM -------------------------------------------------------

MJM_opt <- function(x, y, start = NULL, eps = 0.0001, maxit = 1500,
                    nu = c("lambda" = 0.1, "gamma" = 0.1, "mu" = 1,
                           "sigma" = 1, "alpha" = 1),
                    opt_long = TRUE, alpha.eps = 0.001, par_trace = FALSE,
                    verbose = FALSE, update_nu = FALSE, update_tau = FALSE,
                    nocheck_logpost = FALSE, iwls_sigma = TRUE, coll = FALSE,
                    std_surv = TRUE, ...) {

  if(!is.null(start))
    x <- bamlss::set.starting.values(x, start)
  nsubj <- attr(y, "nsubj")
  gq_weights <- attr(y, "gq_weights")
  n_w <- length(gq_weights)
  take_last <- attr(y, "take_last")
  take_last_l <- attr(y, "take_last_l")
  status <- attr(y, "status")
  survtime <- y[[1]][, "time"][take_last]
  nmarker <- attr(y, "nmarker")
  marker <- attr(y, "marker")
  logLik <- NULL
  edf <- 0
  get.all.par <- utils::getFromNamespace("get.all.par", "bamlss")
  get.hessian <- utils::getFromNamespace("get.hessian", "bamlss")

  if(length(nu) != 5) {
    stop("Please provide the optimizing step length for each predictor.")
  }
  if(is.null(names(nu)))
    names(nu) <- c("lambda", "gamma", "mu", "sigma", "alpha")

  ## Set alpha/mu/sigma intercept starting value.
  eta_timegrid_alpha <- 0
  if(length(x$alpha$smooth.construct)) {
    for(j in names(x$alpha$smooth.construct)) {
      if (j == "model.matrix" & is.null(start)) {
        x$alpha$smooth.construct[[j]]$state$parameters[seq_len(nmarker)] <-
          .Machine$double.eps
        x$alpha$smooth.construct[[j]]$state$fitted.values <-
          x$alpha$smooth.construct[[j]]$X %*%
          x$alpha$smooth.construct[[j]]$state$parameters
        x$alpha$smooth.construct[[j]]$nu <- nu["alpha"]
      }
      b <- get.par(x$alpha$smooth.construct[[j]]$state$parameters, "b")
      eta_grid_sj <- drop(x$alpha$smooth.construct[[j]]$Xgrid %*% b)
      x$alpha$smooth.construct[[j]]$state$fitted_timegrid <- eta_grid_sj
      eta_timegrid_alpha <- eta_timegrid_alpha + eta_grid_sj
      edf <- edf + x$alpha$smooth.construct[[j]]$state$edf
    }
  }
  eta_timegrid_mu <- 0
  eta_T_mu <- 0
  if(length(x$mu$smooth.construct)) {
    for(j in names(x$mu$smooth.construct)) {
      if (j == "model.matrix" & is.null(start)) {
        x$mu$smooth.construct[[j]]$state$parameters[seq_len(nmarker)] <-
          tapply(y[[1]][, "obs"], marker, mean, na.rm = TRUE)
        x$mu$smooth.construct[[j]]$state$fitted.values <-
          x$mu$smooth.construct[[j]]$X %*%
          x$mu$smooth.construct[[j]]$state$parameters
      }
      b <- get.par(x$mu$smooth.construct[[j]]$state$parameters, "b")
      if(inherits(x$mu$smooth.construct[[j]], "pcre2.random.effect")){
        eta_grid_sj <- rep(0, nrow(x$mu$smooth.construct[[j]]$Xgrid))
        eta_T_sj <- rep(0, nrow(x$mu$smooth.construct[[j]]$XT))
      } else {
        eta_grid_sj <- drop(x$mu$smooth.construct[[j]]$Xgrid %*% b)
        eta_T_sj <- drop(x$mu$smooth.construct[[j]]$XT %*% b)
      }
      x$mu$smooth.construct[[j]]$state$fitted_timegrid <- eta_grid_sj
      x$mu$smooth.construct[[j]]$state$fitted_T <- eta_T_sj
      eta_timegrid_mu <- eta_timegrid_mu + eta_grid_sj
      eta_T_mu <- eta_T_mu + eta_T_sj
      x$mu$smooth.construct[[j]]$nu <- nu["mu"]
      edf <- edf + x$mu$smooth.construct[[j]]$state$edf
    }
  }

  eta_timegrid_lambda <- 0
  if(length(x$lambda$smooth.construct)) {
    for(j in names(x$lambda$smooth.construct)) {
      b <- get.par(x$lambda$smooth.construct[[j]]$state$parameters, "b")
      eta_sj <- drop(x$lambda$smooth.construct[[j]]$Xgrid %*% b)
      x$lambda$smooth.construct[[j]]$state$fitted_timegrid <- eta_sj
      x$lambda$smooth.construct[[j]]$nu <- nu[["lambda"]]
      eta_timegrid_lambda <- eta_timegrid_lambda + eta_sj
      edf <- edf + x$lambda$smooth.construct[[j]]$state$edf
    }
  }

  if(!is.null(x$sigma$smooth.construct$model.matrix)) {
    if(is.null(start)) {
      x$sigma$smooth.construct$model.matrix$state$parameters[
        seq_len(nmarker)] <- tapply(y[[1]][, "obs"], marker,
                                    function(x) log(sd(x, na.rm = TRUE)))
    }
    x$sigma$smooth.construct$model.matrix$state$fitted.values <-
      x$sigma$smooth.construct$model.matrix$X %*%
      x$sigma$smooth.construct$model.matrix$state$parameters
    x$sigma$smooth.construct$model.matrix$nu <- nu[["sigma"]]
    edf <- edf + x$sigma$smooth.construct$model.matrix$state$edf
  }

  if (length(x$gamma$smooth.construct)) {
    for (j in names(x$gamma$smooth.construct)) {
      edf <- edf + x$gamma$smooth.construct[[j]]$state$edf
      x$gamma$smooth.construct[[j]]$nu <- nu[["gamma"]]
    }
  }

  get.eta <- utils::getFromNamespace("get.eta", "bamlss")
  eta <- get.eta(x, expand = FALSE)

  # Standardizing the Survival design matrices
  eta_timegrid_mu_unst <- eta_timegrid_mu
  eta_T_mu_unst <- eta_T_mu
  if (!is.null(start) && std_surv) {
    x_std <- matrix(eta_timegrid_mu, ncol = nmarker)
    long_bar <- colMeans(x_std)
    long_sds <- apply(x_std, 2, sd)
  } else {
    long_bar <- rep(0, nmarker)
    long_sds <- rep(1, nmarker)
  }
  attr(eta, "std_long") <- list("long_bar" = long_bar,
                                "long_sds" = long_sds)
  eta_timegrid_mu <-
    (eta_timegrid_mu_unst - rep(long_bar, each = nsubj * n_w)) /
    rep(long_sds, each = nsubj * n_w)
  eta_T_mu <- (eta_T_mu_unst - rep(long_bar, each = nsubj)) /
    rep(long_sds, each = nsubj)

  eta_timegrid_long <- rowSums(matrix(eta_timegrid_alpha*eta_timegrid_mu,
                                      nrow = nsubj*n_w, ncol = nmarker))
  eta_timegrid <- eta_timegrid_lambda + eta_timegrid_long


  ## Extract current value of the log-posterior.
  # Define here so data is available in the function environment
  # Used for updating nu / verbose
  get_LogLik <- function(eta_timegrid, eta_T_mu, eta) {


    eta_T_long <- rowSums(matrix(eta$alpha*eta_T_mu, nrow = nsubj,
                                 ncol = nmarker))
    eta_T <- eta$lambda + eta$gamma + eta_T_long
    sum_Lambda <- drop(
      crossprod(survtime/2 * exp(eta$gamma),
                colSums(gq_weights*matrix(exp(eta_timegrid),
                                          ncol = nsubj,
                                          nrow = length(gq_weights)))))
    logLik <- drop(crossprod(status, eta_T)) - sum_Lambda +
      sum(dnorm(y[[1]][, "obs"], mean = eta$mu, sd = exp(eta$sigma),
                log = TRUE))

    return(logLik)

  }

  # Compare if update increases logPost
  get.log.prior <- utils::getFromNamespace("get.log.prior", "bamlss")
  LogPrioOLD <- get.log.prior(x)
  LogLikOLD <- get_LogLik(eta_timegrid, eta_T_mu, eta)
  etaUP <- eta

  # For algorithm
  eps0 <- eps0_surv <- eps0_long <- eps + 1
  eta0_surv <- do.call("cbind", eta[c("lambda", "gamma")])
  eta0_alpha <- matrix(eta$alpha, nrow = nsubj, ncol = nmarker)
  eta0 <- cbind(eta0_surv, eta0_alpha)
  eta0_long <- do.call("cbind", eta[c("mu", "sigma")])
  iter <- 0
  alpha_update <- if(is.null(start)) FALSE else TRUE

  # Parameter trace to follow the estimated parameter values
  if (par_trace) {
    it_param <- list()
  }

  # Collinearity measures for the alpha predictors
  if (coll) {
    it_coll <- list()
  }

  # Updating the predictors
  while((eps0 > eps) & (iter < maxit)) {
    ## (1) update lambda.
    for(j in names(x$lambda$smooth.construct)) {
      state <- update_mjm_lambda(x$lambda$smooth.construct[[j]], y = y,
                                 eta = eta, eta_timegrid = eta_timegrid,
                                 eta_T_mu = eta_T_mu,
                                 survtime = survtime, update_nu = update_nu,
                                 get_LogLik = get_LogLik,
                                 update_tau = update_tau, edf = edf, ...)
      eta_timegrid_lambdaUP <- eta_timegrid_lambda -
        x$lambda$smooth.construct[[j]]$state$fitted_timegrid +
        state$fitted_timegrid
      eta_timegridUP <- eta_timegrid_lambdaUP +
        eta_timegrid_long
      etaUP$lambda <- eta$lambda - fitted(x$lambda$smooth.construct[[j]]$state) +
        fitted(state)

      # Check for update
      LogPrioUP <- LogPrioOLD -
        x$lambda$smooth.construct[[j]]$prior(
          x$lambda$smooth.construct[[j]]$state$parameters) +
        x$lambda$smooth.construct[[j]]$prior(state$parameters)
      LogLikUP <- get_LogLik(eta_timegridUP, eta_T_mu, etaUP)
      if (nocheck_logpost || LogPrioUP + LogLikUP > LogPrioOLD + LogLikOLD ||
          iter == 0) {
        LogPrioOLD <- LogPrioUP
        LogLikOLD <- LogLikUP
        eta_timegrid_lambda <- eta_timegrid_lambdaUP
        eta_timegrid <- eta_timegridUP
        eta$lambda <- etaUP$lambda
        x$lambda$smooth.construct[[j]]$state <- state
      } else {
        etaUP$lambda <- eta$lambda
      }
    }

    ## (2) update gamma.
    if(length(x$gamma$smooth.construct)) {
      for(j in seq_along(x$gamma$smooth.construct)) {
        state <- update_mjm_gamma(x$gamma$smooth.construct[[j]], y = y,
                                  eta = eta, eta_timegrid = eta_timegrid,
                                  eta_T_mu = eta_T_mu,
                                  survtime = survtime, update_nu = update_nu,
                                  get_LogLik = get_LogLik,
                                  update_tau = update_tau, edf = edf, ...)
        etaUP$gamma <- eta$gamma - fitted(x$gamma$smooth.construct[[j]]$state) +
          fitted(state)

        # Check for update
        LogPrioUP <- LogPrioOLD -
          x$gamma$smooth.construct[[j]]$prior(
            x$gamma$smooth.construct[[j]]$state$parameters) +
          x$gamma$smooth.construct[[j]]$prior(state$parameters)
        LogLikUP <- get_LogLik(eta_timegrid, eta_T_mu, etaUP)
        if (nocheck_logpost || LogPrioUP + LogLikUP > LogPrioOLD + LogLikOLD ||
            iter == 0) {
          LogPrioOLD <- LogPrioUP
          LogLikOLD <- LogLikUP
          eta$gamma <- etaUP$gamma
          x$gamma$smooth.construct[[j]]$state <- state
        } else {
          etaUP$gamma <- eta$gamma
        }
      }
    }

    if (!alpha_update && max(c(eps0_surv, eps0_long)) < alpha.eps) {
      alpha_update <- 1
      if (verbose) {
        cat("It ", iter, "-- Start Alpha Update", "\n")
      }
    }
    if (opt_long) {
      ## (3) update alpha.
      if(alpha_update) {
        if(alpha_update == 1 && std_surv) {
          x_std <- matrix(eta_timegrid_mu, ncol = nmarker)
          long_bar <- colMeans(x_std)
          long_sds <- apply(x_std, 2, sd)
          attr(eta, "std_long") <- list(
            "long_bar" = long_bar,
            "long_sds" = long_sds
          )
          if(is.null(start)) {
            eta_timegrid_mu <-
              (eta_timegrid_mu_unst - rep(long_bar, each = nsubj * n_w)) /
              rep(long_sds, each = nsubj * n_w)
            eta_T_mu <- (eta_T_mu_unst - rep(long_bar, each = nsubj)) /
              rep(long_sds, each = nsubj)
          }
        }
        if(length(x$alpha$smooth.construct)) {
          for(j in seq_along(x$alpha$smooth.construct)) {
            state <- update_mjm_alpha(x$alpha$smooth.construct[[j]], y = y,
                                      eta = eta, eta_timegrid = eta_timegrid,
                                      eta_timegrid_lambda = eta_timegrid_lambda,
                                      eta_timegrid_mu = eta_timegrid_mu,
                                      eta_timegrid_alpha = eta_timegrid_alpha,
                                      eta_T_mu = eta_T_mu, survtime = survtime,
                                      update_nu = update_nu,
                                      get_LogLik = get_LogLik,
                                      update_tau = update_tau, edf = edf,
                                      coll = coll, ...)
            etaUP$alpha <- eta$alpha -
              drop(fitted(x$alpha$smooth.construct[[j]]$state)) +
              fitted(state)
            eta_timegrid_alphaUP <- eta_timegrid_alpha -
              x$alpha$smooth.construct[[j]]$state$fitted_timegrid +
              state$fitted_timegrid
            eta_timegrid_longUP <- rowSums(matrix(eta_timegrid_alphaUP *
                                                  eta_timegrid_mu,
                                                nrow = nsubj*n_w,
                                                ncol = nmarker))
            eta_timegridUP <- eta_timegrid_lambda + eta_timegrid_longUP

            # Check for update
            LogPrioUP <- LogPrioOLD -
              x$alpha$smooth.construct[[j]]$prior(
                x$alpha$smooth.construct[[j]]$state$parameters) +
              x$alpha$smooth.construct[[j]]$prior(state$parameters)
            LogLikUP <- get_LogLik(eta_timegridUP, eta_T_mu, etaUP)
            if (nocheck_logpost ||
                LogPrioUP + LogLikUP > LogPrioOLD + LogLikOLD ||
                alpha_update == 1) {
              LogPrioOLD <- LogPrioUP
              LogLikOLD <- LogLikUP
              eta_timegrid_alpha <- eta_timegrid_alphaUP
              eta_timegrid_long <- eta_timegrid_longUP
              eta_timegrid <- eta_timegridUP
              eta$alpha <- etaUP$alpha
              x$alpha$smooth.construct[[j]]$state <- state
              alpha_update <- alpha_update + 1
            } else {
              etaUP$alpha <- eta$alpha
            }

            # Collinearity measures for alpha
            if (coll) {
              it_coll[[iter]] <- state$coll
            }
          }
        }
      }

      ## (4) update mu.
      if(length(x$mu$smooth.construct)) {
        for(j in seq_along(x$mu$smooth.construct)) {
          state <- update_mjm_mu(x$mu$smooth.construct[[j]], y = y,
                                 eta = eta, eta_timegrid = eta_timegrid,
                                 eta_timegrid_lambda = eta_timegrid_lambda,
                                 eta_timegrid_mu_unst = eta_timegrid_mu_unst,
                                 eta_timegrid_alpha = eta_timegrid_alpha,
                                 eta_T_mu_unst = eta_T_mu_unst,
                                 survtime = survtime,
                                 get_LogLik = get_LogLik, update_nu = update_nu,
                                 update_tau = update_tau, edf = edf)
          etaUP$mu <- eta$mu -
            drop(fitted(x$mu$smooth.construct[[j]]$state)) + fitted(state)
          eta_timegrid_mu_unstUP <- eta_timegrid_mu_unst -
            x$mu$smooth.construct[[j]]$state$fitted_timegrid +
            state$fitted_timegrid
          eta_timegrid_muUP <-
            (eta_timegrid_mu_unstUP - rep(long_bar, each = nsubj * n_w)) /
            rep(long_sds, each = nsubj * n_w)
          eta_timegrid_longUP <- rowSums(matrix(eta_timegrid_alpha *
                                                eta_timegrid_muUP,
                                              nrow = nsubj*n_w, ncol = nmarker))
          eta_timegridUP <- eta_timegrid_lambda + eta_timegrid_longUP
          eta_T_mu_unstUP <- eta_T_mu_unst -
            x$mu$smooth.construct[[j]]$state$fitted_T +
            state$fitted_T
          eta_T_muUP <-  (eta_T_mu_unst - rep(long_bar, each = nsubj)) /
            rep(long_sds, each = nsubj)

          # Check for update
          LogPrioUP <- LogPrioOLD -
            x$mu$smooth.construct[[j]]$prior(
              x$mu$smooth.construct[[j]]$state$parameters) +
            x$mu$smooth.construct[[j]]$prior(state$parameters)
          LogLikUP <- get_LogLik(eta_timegridUP, eta_T_muUP, etaUP)
          if (nocheck_logpost ||
              LogPrioUP + LogLikUP > LogPrioOLD + LogLikOLD || iter == 0) {
            LogPrioOLD <- LogPrioUP
            LogLikOLD <- LogLikUP
            eta_timegrid_mu_unst <- eta_timegrid_mu_unstUP
            eta_timegrid_mu <- eta_timegrid_muUP
            eta_timegrid_long <- eta_timegrid_longUP
            eta_timegrid <- eta_timegridUP
            eta_T_mu_unst <- eta_T_mu_unstUP
            eta_T_mu <- eta_T_muUP
            eta$mu <- etaUP$mu
            x$mu$smooth.construct[[j]]$state <- state
          } else {
            etaUP$mu <- eta$mu
          }
        }
      }

      ## (5) update sigma.
      if(length(x$sigma$smooth.construct)) {
        for(j in seq_along(x$sigma$smooth.construct)) {
          state <- update_mjm_sigma(x$sigma$smooth.construct[[j]], y = y,
                                    eta = eta,
                                    eta_timegrid = eta_timegrid,
                                    eta_T_mu = eta_T_mu,
                                    get_LogLik = get_LogLik,
                                    survtime = survtime, update_nu = update_nu,
                                    update_tau = update_tau, edf = edf,
                                    iwls = iwls_sigma, ...)
          etaUP$sigma <- eta$sigma -
            drop(fitted(x$sigma$smooth.construct[[j]]$state)) +
            fitted(state)

          # Check for update
          LogPrioUP <- LogPrioOLD -
            x$sigma$smooth.construct[[j]]$prior(
              x$sigma$smooth.construct[[j]]$state$parameters) +
            x$sigma$smooth.construct[[j]]$prior(state$parameters)
          LogLikUP <- get_LogLik(eta_timegrid, eta_T_mu, etaUP)
          if (nocheck_logpost ||
              LogPrioUP + LogLikUP > LogPrioOLD + LogLikOLD || iter == 0) {
            LogPrioOLD <- LogPrioUP
            LogLikOLD <- LogLikUP
            eta$sigma <- etaUP$sigma
            x$sigma$smooth.construct[[j]]$state <- state
          } else {
            etaUP$sigma <- eta$sigma
          }
        }
      }
    }


    if (verbose) {
      # Likelihood calculation

      logLik <- get_LogLik(eta_timegrid = eta_timegrid,
                           eta_T_mu = eta_T_mu,
                           eta = eta)

      cat("It ", iter,", LogLik ", logLik, ", Post",
          as.numeric(logLik + get.log.prior(x)), "\n")
    }


    # Stopping criterion
    iter <- iter + 1
    eta1_surv <- do.call("cbind", eta[c("lambda", "gamma")])
    eta1_alpha <- matrix(eta$alpha, nrow = nsubj, ncol = nmarker)
    eta1_long <- do.call("cbind", eta[c("mu", "sigma")])
    eps0_surv <- mean(abs((eta1_surv - eta0_surv) / eta1_surv), na.rm = TRUE)
    eps0_alpha <- if (alpha_update) {
      mean(abs((eta1_alpha - eta0_alpha) / eta1_alpha), na.rm = TRUE)
    } else {
      eps + 1
    }
    eps0_long <- mean(abs((eta1_long - eta0_long) / eta1_long), na.rm = TRUE)
    eps0 <- mean(c(eps0_surv, eps0_alpha, eps0_long))

    eta0_surv <- eta1_surv
    eta0_alpha <- eta1_alpha
    eta0_long <- eta1_long

    # Parameter trace to follow the estimated parameter values
    if (par_trace) {
      it_param[[iter]] <- get.all.par(x)
    }

  }

  # Log-Posterior ausrechnen und ausgeben
  if (is.null(logLik)) {
    logLik <- get_LogLik(eta_timegrid = eta_timegrid,
                         eta_T_mu = eta_T_mu,
                         eta = eta)
  }
  logPost <- as.numeric(logLik + get.log.prior(x))
  return(list("fitted.values" = eta, "parameters" = get.all.par(x),
              "logLik" = logLik, "logPost" = logPost,
              "hessian" = get.hessian(x),
              "converged" = iter < maxit,
              "scaling" = attr(eta, "std_long"),
              "par_trace" = if (par_trace) it_param else NULL,
              "coll" = if (coll) it_coll else NULL))

}
