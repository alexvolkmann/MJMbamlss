
# Updating lambda predictor -----------------------------------------------


update_mjm_lambda <- function(x, y, eta, eta_timegrid, eta_T_mu, survtime,
                              update_nu, get_LogLik, update_tau, edf, ...)
{

  b <- bamlss::get.state(x, "b")
  b_p <- length(b)
  nu <- x$nu
  matrix_inv <- utils::getFromNamespace("matrix_inv", "bamlss")
  sum_diag <- utils::getFromNamespace("sum_diag", "bamlss")
  get.ic2 <- utils::getFromNamespace("get.ic2", "bamlss")
  tau2.optim <- utils::getFromNamespace("tau2.optim", "bamlss")

  int_i <- survint_C(pred = "lambda", pre_fac = exp(eta$gamma),
                      omega = exp(eta_timegrid),
                      int_vec = x$Xgrid, weights = attr(y, "gq_weights"),
                      survtime = survtime)


  # Status from MJM_transform
  # XT from sm_time_transform_mjm()
  x_score0 <- drop(attr(y, "status") %*% x$XT) - int_i$score_int
  x_H0 <- matrix(int_i$hess_int, ncol = b_p)

  # Updating Taus
  if(update_tau && (x$state$do.optim || !x$fixed)) {
    env <- new.env()
    par <- x$state$parameters
    edf0 <- edf - x$state$edf
    status <- attr(y, "status")

    # Function for Updating Taus
    # Defined here so that it accesses x_score and x_H
    tau_update_opt <- function(tau2) {
      par[x$pid$tau2] <- tau2
      x_score <- x_score0 + x$grad(score = NULL, par, full = FALSE)
      x_H <- x_H0 + x$hess(score = NULL, par, full = FALSE)
      Sigma <- matrix_inv(x_H, index = NULL)
      Hs <- Sigma %*% x_score
      if(update_nu) {
        # Function for updating step length
        # Defined here so that it accesses Hs and par and fitted from parents
        nu_update_opt <- function(nu) {
          b2 <- drop(b + nu * Hs)
          par[x$pid$b] <- b2
          fit <- drop(x$X %*% b2)
          fitted_timegrid <- drop(x$Xgrid %*% b2)
          eta$lambda <- eta$lambda - fitted(x$state) + fit
          eta_timegrid <- eta_timegrid - x$state$fitted_timegrid +
            fitted_timegrid
          LogPost <- get_LogLik(eta_timegrid, eta_T_mu, eta) + x$prior(par)
          return(-1 * LogPost)
        }
        nu <- optimize(f = nu_update_opt, interval = c(0, 1))$minimum
      }

      # Calculate the Information Criterion with given edf and nu
      b2 <- drop(b + nu * Hs)
      fit <- drop(x$X %*% b2)
      fitted_timegrid <- drop(x$Xgrid %*% b2)
      eta$lambda <- eta$lambda - fitted(x$state) + fit
      eta_timegrid <- eta_timegrid - x$state$fitted_timegrid +
        fitted_timegrid
      edf1 <- sum_diag(x_H0 %*% Sigma)
      edf <- edf0 + edf1
      logLik <- get_LogLik(eta_timegrid, eta_T_mu, eta)
      ic <- get.ic2(logLik, edf, n = length(eta$mu), type = "AICc")

      # First time running the function or new minimum
      if (is.null(env$ic_val) || (ic < env$ic_val)) {

        # Prepare the current state for output
        x$state$parameters[x$pid$b] <- b2
        x$state$parameters[x$pid$tau2] <- tau2
        x$state$fitted_timegrid <- fitted_timegrid
        x$state$fitted.values <- fit
        x$state$edf <- edf1
        x$state$hessian <- x_H

        assign("ic_val", ic, envir = env)
        assign("state", x$state, envir = env)
      }

      # For optimizer
      return(ic)
    }

    ic0 <- tau_update_opt(get.state(x, "tau2"))
    tau2 <- tau2.optim(tau_update_opt, start = get.state(x, "tau2"))
    return(env$state)

  }

  ## Newton-Raphson.
  x_score <- x_score0 + x$grad(score = NULL, x$state$parameters, full = FALSE)
  x_H <- x_H0 + x$hess(score = NULL, x$state$parameters, full = FALSE)
  Sigma <- matrix_inv(x_H, index = NULL)
  Hs <- Sigma %*% x_score

  if(update_nu) {
    par <- x$state$parameters
    # Function for updating step length
    # Defined here so that it accesses Hs and par and fitted from parents
    nu_update_opt <- function(nu) {
      b2 <- drop(b + nu * Hs)
      par[x$pid$b] <- b2
      fit <- drop(x$X %*% b2)
      fitted_timegrid <- drop(x$Xgrid %*% b2)
      eta$lambda <- eta$lambda - fitted(x$state) + fit
      eta_timegrid <- eta_timegrid - x$state$fitted_timegrid +
        fitted_timegrid
      LogPost <- get_LogLik(eta_timegrid, eta_T_mu, eta) + x$prior(par)
      return(-1 * LogPost)
    }
    nu <- optimize(f = nu_update_opt, interval = c(0, 1))$minimum
  }
  b <- b + nu * Hs

  x$state$parameters[x$pid$b] <- b
  x$state$fitted_timegrid <- drop(x$Xgrid %*% b)
  x$state$fitted.values <- drop(x$X %*% b)
  x$state$hessian <- x_H
  x$state$edf <- sum_diag(x_H0 %*% Sigma)
  return(x$state)

}



# Updating gamma predictor ------------------------------------------------


update_mjm_gamma <- function(x, y, eta, eta_timegrid, eta_T_mu, survtime,
                             update_nu, get_LogLik, update_tau, edf, ...) {

  b <- bamlss::get.state(x, "b")
  b_p <- length(b)
  take_last <- attr(y, "take_last")
  exp_eta_gamma <- exp(eta$gamma)
  nu <- x$nu
  matrix_inv <- utils::getFromNamespace("matrix_inv", "bamlss")
  sum_diag <- utils::getFromNamespace("sum_diag", "bamlss")
  get.ic2 <- utils::getFromNamespace("get.ic2", "bamlss")
  tau2.optim <- utils::getFromNamespace("tau2.optim", "bamlss")

  int_i <- survint_C(pred = "gamma", pre_fac = exp_eta_gamma, pre_vec = x$X,
                      omega = exp(eta_timegrid),
                      weights = attr(y, "gq_weights"),
                      survtime = survtime)
  x_score0 <- drop(attr(y, "status") %*% x$X) - int_i$score_int
  x_H0 <- matrix(int_i$hess_int, ncol = b_p)

  # Updating Taus
  if(update_tau && (x$state$do.optim || !x$fixed)) {
    env <- new.env()
    par <- x$state$parameters
    edf0 <- edf - x$state$edf
    status <- attr(y, "status")

    # Function for Updating Taus
    # Defined here so that it accesses x_score and x_H
    tau_update_opt <- function(tau2) {
      par[x$pid$tau2] <- tau2
      x_score <- x_score0 + x$grad(score = NULL, par, full = FALSE)
      x_H <- x_H0 + x$hess(score = NULL, par, full = FALSE)
      Sigma <- matrix_inv(x_H, index = NULL)
      Hs <- Sigma %*% x_score
      if(update_nu) {
        # Function for updating step length
        # Defined here so that it accesses Hs and par and fitted from parents
        nu_update_opt <- function(nu) {
          b2 <- drop(b + nu * Hs)
          par[x$pid$b] <- b2
          fit <- drop(x$X %*% b2)
          eta$gamma <- eta$gamma - fitted(x$state) + fit
          LogPost <- get_LogLik(eta_timegrid, eta_T_mu, eta) + x$prior(par)
          return(-1 * LogPost)
        }
        nu <- optimize(f = nu_update_opt, interval = c(0, 1))$minimum
      }

      # Calculate the Information Criterion with given edf and nu
      b2 <- drop(b + nu * Hs)
      fit <- drop(x$X %*% b2)
      eta$gamma <- eta$gamma - fitted(x$state) + fit
      edf1 <- sum_diag(x_H0 %*% Sigma)
      edf <- edf0 + edf1
      logLik <- get_LogLik(eta_timegrid, eta_T_mu, eta)
      ic <- get.ic2(logLik, edf, n = length(eta$mu), type = "AICc")

      # First time running the function or new minimum
      if (is.null(env$ic_val) || (ic < env$ic_val)) {

        # Prepare the current state for output
        x$state$parameters[x$pid$b] <- b2
        x$state$parameters[x$pid$tau2] <- tau2
        x$state$fitted.values <- fit
        x$state$edf <- edf1
        x$state$hessian <- x_H

        assign("ic_val", ic, envir = env)
        assign("state", x$state, envir = env)
      }

      # For optimizer
      return(ic)
    }

    ic0 <- tau_update_opt(get.state(x, "tau2"))
    tau2 <- tau2.optim(tau_update_opt, start = get.state(x, "tau2"))
    return(env$state)

  }

  x_score <- x_score0 + x$grad(score = NULL, x$state$parameters, full = FALSE)
  x_H <- x_H0 + x$hess(score = NULL, x$state$parameters, full = FALSE)
  Sigma <- matrix_inv(x_H, index = NULL)
  Hs <- Sigma %*% x_score

  if(update_nu) {
    par <- x$state$parameters
    # Function for updating step length
    # Defined here so that it accesses Hs and par and fitted from parents
    nu_update_opt <- function(nu) {
      b2 <- drop(b + nu * Hs)
      par[x$pid$b] <- b2
      fit <- drop(x$X %*% b2)
      eta$gamma <- eta$gamma - fitted(x$state) + fit
      LogPost <- get_LogLik(eta_timegrid, eta_T_mu, eta) + x$prior(par)
      return(-1 * LogPost)
    }
    nu <- optimize(f = nu_update_opt, interval = c(0, 1))$minimum
  }
  b <- b + nu * Hs

  x$state$parameters[seq_len(b_p)] <- b
  x$state$fitted.values <- drop(x$X %*% b)
  x$state$hessian <- x_H
  x$state$edf <- sum_diag(x_H0 %*% Sigma)
  return(x$state)

}


# Updating alpha predictor ------------------------------------------------


update_mjm_alpha <- function(x, y, eta, eta_timegrid, eta_timegrid_lambda,
                             eta_timegrid_alpha, eta_timegrid_mu,
                             eta_T_mu, survtime, update_nu, get_LogLik,
                             update_tau, edf, coll, ...) {

  b <- bamlss::get.state(x, "b")
  b_p <- length(b)
  nmarker <- attr(y, "nmarker")
  status <- attr(y, "status")
  nsubj <- attr(y, "nsubj")
  n_w <- length(attr(y, "gq_weights"))
  nu <- x$nu
  matrix_inv <- utils::getFromNamespace("matrix_inv", "bamlss")
  sum_diag <- utils::getFromNamespace("sum_diag", "bamlss")
  get.ic2 <- utils::getFromNamespace("get.ic2", "bamlss")
  tau2.optim <- utils::getFromNamespace("tau2.optim", "bamlss")

  int_i <- survint_C(pred = "long", pre_fac = exp(eta$gamma),
                      omega = exp(eta_timegrid),
                      int_fac = eta_timegrid_mu, int_vec = x$Xgrid,
                      weights = attr(y, "gq_weights"),
                      survtime = survtime)
  delta <- rep(attr(y, "status"), nmarker)
  x_score0 <- drop(t(delta * x$XT) %*% eta_T_mu) - int_i$score_int
  x_H0 <- matrix(int_i$hess_int, ncol = b_p)

  # Updating Taus
  if(update_tau && (x$state$do.optim || !x$fixed)) {
    env <- new.env()
    par <- x$state$parameters
    edf0 <- edf - x$state$edf

    # Function for Updating Taus
    # Defined here so that it accesses x_score and x_H
    tau_update_opt <- function(tau2) {
      par[x$pid$tau2] <- tau2
      x_score <- x_score0 + x$grad(score = NULL, par, full = FALSE)
      x_H <- x_H0 + x$hess(score = NULL, par, full = FALSE)
      Sigma <- matrix_inv(x_H, index = NULL)
      Hs <- Sigma %*% x_score
      if(update_nu) {
        # Function for updating step length
        # Defined here so that it accesses Hs and par and fitted from parents
        nu_update_opt <- function(nu) {
          b2 <- drop(b + nu * Hs)
          par[x$pid$b] <- b2
          fit <- drop(x$X %*% b2)
          fitted_timegrid <- drop(x$Xgrid %*% b2)
          eta$alpha <- eta$alpha - fitted(x$state) + fit
          eta_timegrid_alpha <- eta_timegrid_alpha -
            x$state$fitted_timegrid +fitted_timegrid
          eta_timegrid_long <- rowSums(matrix(eta_timegrid_alpha *
                                                eta_timegrid_mu,
                                              nrow = nsubj*n_w,
                                              ncol = nmarker))
          eta_timegrid <- eta_timegrid_lambda + eta_timegrid_long
          LogPost <- get_LogLik(eta_timegrid, eta_T_mu, eta) + x$prior(par)
          return(-1 * LogPost)
        }
        nu <- optimize(f = nu_update_opt, interval = c(0, 1))$minimum
      }

      # Calculate the Information Criterion with given edf and nu
      b2 <- drop(b + nu * Hs)
      fit <- drop(x$X %*% b2)
      fitted_timegrid <- drop(x$Xgrid %*% b2)
      eta$alpha <- eta$alpha - fitted(x$state) + fit
      eta_timegrid_alpha <- eta_timegrid_alpha -
        x$state$fitted_timegrid +fitted_timegrid
      eta_timegrid_long <- rowSums(matrix(eta_timegrid_alpha *
                                            eta_timegrid_mu,
                                          nrow = nsubj*n_w,
                                          ncol = nmarker))
      eta_timegrid <- eta_timegrid_lambda + eta_timegrid_long
      edf1 <- sum_diag(x_H0 %*% Sigma)
      edf <- edf0 + edf1
      logLik <- get_LogLik(eta_timegrid, eta_T_mu, eta)
      ic <- get.ic2(logLik, edf, n = length(eta$mu), type = "AICc")

      # First time running the function or new minimum
      if (is.null(env$ic_val) || (ic < env$ic_val)) {

        # Prepare the current state for output
        x$state$parameters[x$pid$b] <- b2
        x$state$parameters[x$pid$tau2] <- tau2
        x$state$fitted_timegrid <- fitted_timegrid
        x$state$fitted.values <- fit
        x$state$edf <- edf1
        x$state$hessian <- x_H

        assign("ic_val", ic, envir = env)
        assign("state", x$state, envir = env)
      }

      # For optimizer
      return(ic)
    }

    ic0 <- tau_update_opt(get.state(x, "tau2"))
    tau2 <- tau2.optim(tau_update_opt, start = get.state(x, "tau2"))
    return(env$state)

  }


  # Newton Raphson
  x_score <- x_score0 + x$grad(score = NULL, x$state$parameters, full = FALSE)
  p_H <- x$hess(score = NULL, x$state$parameters, full = FALSE)
  x_H <- x_H0 + p_H
  Sigma <- matrix_inv(x_H, index = NULL)
  Hs <- Sigma %*% x_score

  if(update_nu) {
    par <- x$state$parameters
    # Function for updating step length
    # Defined here so that it accesses Hs and par and fitted from parents
    nu_update_opt <- function(nu) {
      b2 <- drop(b + nu * Hs)
      par[x$pid$b] <- b2
      fit <- drop(x$X %*% b2)
      fitted_timegrid <- drop(x$Xgrid %*% b2)
      eta$alpha <- eta$alpha - fitted(x$state) + fit
      eta_timegrid_alpha <- eta_timegrid_alpha -
        x$state$fitted_timegrid +fitted_timegrid
      eta_timegrid_long <- rowSums(matrix(eta_timegrid_alpha *
                                            eta_timegrid_mu,
                                          nrow = nsubj*n_w,
                                          ncol = nmarker))
      eta_timegrid <- eta_timegrid_lambda + eta_timegrid_long
      LogPost <- get_LogLik(eta_timegrid, eta_T_mu, eta) + x$prior(par)
      return(-1 * LogPost)
    }
    nu <- optimize(f = nu_update_opt, interval = c(0, 1))$minimum
  }
  b <- b + nu * Hs

  x$state$parameters[seq_len(b_p)] <- b
  x$state$fitted_timegrid <- drop(x$Xgrid %*% b)
  x$state$fitted.values <- drop(x$X %*% b)
  x$state$hessian <- x_H
  x$state$edf <- sum_diag(x_H0 %*% Sigma)
  if(coll) {
    x$state$coll <- list("X" = eta_T_mu, "I" = int_i$hess_int, "P" = p_H)
  }

  return(x$state)

}



# Updating mu predictor ---------------------------------------------------


update_mjm_mu <- function(x, y, eta, eta_timegrid, eta_timegrid_lambda,
                          eta_timegrid_alpha, eta_timegrid_mu_unst,
                          eta_T_mu_unst, survtime, update_nu, get_LogLik,
                          update_tau, edf, ...) {

  b <- bamlss::get.state(x, "b")
  b_p <- length(b)
  nmarker <- attr(y, "nmarker")
  delta <- rep(attr(y, "status"), nmarker)
  status <- attr(y, "status")
  nsubj <- attr(y, "nsubj")
  n_w <- length(attr(y, "gq_weights"))
  nu <- x$nu
  long_bar <- attr(eta, "std_long")$long_bar
  long_sds <- attr(eta, "std_long")$long_sds
  matrix_inv <- utils::getFromNamespace("matrix_inv", "bamlss")
  sum_diag <- utils::getFromNamespace("sum_diag", "bamlss")
  get.ic2 <- utils::getFromNamespace("get.ic2", "bamlss")
  tau2.optim <- utils::getFromNamespace("tau2.optim", "bamlss")

  if (any(class(x) == "unc_pcre.random.effect")) {
    int_i <- survint_C(pred = "fpc_re", pre_fac = exp(eta$gamma),
                       omega = exp(eta_timegrid),
                       int_fac = eta_timegrid_alpha /
                         rep(long_sds, each = nsubj*n_w),
                       int_vec = x$Xgrid,
                       weights = attr(y, "gq_weights"),
                       survtime = survtime)
    x_score0 <- drop(
      crossprod(x$X, (y[[1]][, "obs"] - eta$mu) / exp(eta$sigma)^2)  +
        crossprod(delta * x$XT, eta$alpha / rep(long_sds, each = nsubj))) -
      int_i$score_int
    x_H0 <- diag(psi_mat_crossprod(Psi = x, R = 1 / exp(eta$sigma)^2) +
                  int_i$hess_int)
  } else {
    int_i <- survint_C(pred = "long", pre_fac = exp(eta$gamma),
                       omega = exp(eta_timegrid),
                       int_fac = eta_timegrid_alpha /
                         rep(long_sds, each = nsubj*n_w), int_vec = x$Xgrid,
                       weights = attr(y, "gq_weights"),
                       survtime = survtime)
    x_score0 <- drop(
      crossprod(x$X, (y[[1]][, "obs"] - eta$mu) / exp(eta$sigma)^2)  +
        crossprod(delta * x$XT, eta$alpha / rep(long_sds, each = nsubj))) -
      int_i$score_int
    x_H0 <- crossprod(x$X * (1 / exp(eta$sigma)^2), x$X) +
      matrix(int_i$hess_int, ncol = b_p)
  }

  # Updating Taus
  if(update_tau && (x$state$do.optim || !x$fixed)) {
    env <- new.env()
    par <- x$state$parameters
    edf0 <- edf - x$state$edf

    # Function for Updating Taus
    # Defined here so that it accesses x_score and x_H
    tau_update_opt <- function(tau2) {
      par[x$pid$tau2] <- tau2
      x_score <- x_score0 + x$grad(score = NULL, par, full = FALSE)
      x_H <- x_H0 + x$hess(score = NULL, par, full = FALSE)
      Sigma <- matrix_inv(x_H, index = NULL)
      Hs <- Sigma %*% x_score
      if(update_nu) {
        # Function for updating step length
        # Defined here so that it accesses Hs and par and fitted from parents
        nu_update_opt <- function(nu) {
          b2 <- drop(b + nu * Hs)
          if ("random.effect" %in% class(x)) {
            b2 <- b2 - mean(b2)
          }
          par[x$pid$b] <- b2
          fit <- drop(x$X %*% b2)
          fitted_timegrid <- drop(x$Xgrid %*% b2)
          fitted_T <- drop(x$XT %*% b)
          eta$mu <- eta$mu - fitted(x$state) + fit
          eta_timegrid_mu_unst <- eta_timegrid_mu_unst -
            x$state$fitted_timegrid + fitted_timegrid
          eta_timegrid_mu <-
            (eta_timegrid_mu_unst - rep(long_bar, each = nsubj*n_w)) /
            rep(long_sds, each = nsubj*n_w)
          eta_timegrid_long <- rowSums(matrix(eta_timegrid_alpha *
                                                eta_timegrid_mu,
                                              nrow = nsubj*n_w,
                                              ncol = nmarker))
          eta_timegrid <- eta_timegrid_lambda + eta_timegrid_long
          eta_T_mu_unst <- eta_T_mu_unst - x$state$fitted_T + fitted_T
          eta_T_mu <- (eta_T_mu_unst - rep(long_bar, each = nsubj)) /
            rep(long_sds, each = nsubj)
          LogPost <- get_LogLik(eta_timegrid, eta_T_mu, eta) + x$prior(par)
          return(-1 * LogPost)
        }
        nu <- optimize(f = nu_update_opt, interval = c(0, 1))$minimum
      }

      # Calculate the Information Criterion with given edf and nu
      b2 <- drop(b + nu * Hs)
      if ("random.effect" %in% class(x)) {
        b2 <- b2 - mean(b2)
      }
      fit <- drop(x$X %*% b2)
      fitted_timegrid <- drop(x$Xgrid %*% b2)
      fitted_T <- drop(x$XT %*% b)
      eta$mu <- eta$mu - fitted(x$state) + fit
      eta_timegrid_mu_unst <- eta_timegrid_mu_unst -
        x$state$fitted_timegrid + fitted_timegrid
      eta_timegrid_mu <-
        (eta_timegrid_mu_unst - rep(long_bar, each = nsubj*n_w)) /
        rep(long_sds, each = nsubj*n_w)
      eta_timegrid_long <- rowSums(matrix(eta_timegrid_alpha *
                                            eta_timegrid_mu,
                                          nrow = nsubj*n_w,
                                          ncol = nmarker))
      eta_timegrid <- eta_timegrid_lambda + eta_timegrid_long
      eta_T_mu_unst <- eta_T_mu_unst - x$state$fitted_T + fitted_T
      eta_T_mu <- (eta_T_mu_unst - rep(long_bar, each = nsubj)) /
        rep(long_sds, each = nsubj)
      edf1 <- if (any(class(x) == "unc_pcre.random.effect")) {
        sum(diag(x_H0) * diag(Sigma))
      } else {
        sum_diag(x_H0 %*% Sigma)
      }
      edf <- edf0 + edf1
      logLik <- get_LogLik(eta_timegrid, eta_T_mu, eta)
      ic <- get.ic2(logLik, edf, n = length(eta$mu), type = "AICc")

      # First time running the function or new minimum
      if (is.null(env$ic_val) || (ic < env$ic_val)) {

        # Prepare the current state for output
        x$state$parameters[x$pid$b] <- b2
        x$state$parameters[x$pid$tau2] <- tau2
        x$state$fitted_timegrid <- fitted_timegrid
        x$state$fitted.values <- fit
        x$state$fitted_T <- drop(x$XT %*% b)
        x$state$edf <- edf1
        x$state$hessian <- x_H

        assign("ic_val", ic, envir = env)
        assign("state", x$state, envir = env)
      }

      # For optimizer
      return(ic)
    }

    ic0 <- tau_update_opt(get.state(x, "tau2"))
    tau2 <- tau2.optim(tau_update_opt, start = get.state(x, "tau2"))
    return(env$state)

  }


  # Newton Raphson
  x_score <- x_score0 + x$grad(score = NULL, x$state$parameters, full = FALSE)
  x_H <- x_H0 + x$hess(score = NULL, x$state$parameters, full = FALSE)
  Sigma <- matrix_inv(x_H, index = NULL)
  Hs <- Sigma %*% x_score

  if(update_nu) {
    par <- x$state$parameters
    # Function for updating step length
    # Defined here so that it accesses Hs and par and fitted from parents
    nu_update_opt <- function(nu) {
      b2 <- drop(b + nu * Hs)
      if ("random.effect" %in% class(x)) {
        b2 <- b2 - mean(b2)
      }
      par[x$pid$b] <- b2
      fit <- drop(x$X %*% b2)
      fitted_timegrid <- drop(x$Xgrid %*% b2)
      fitted_T <- drop(x$XT %*% b)
      eta$mu <- eta$mu - fitted(x$state) + fit
      eta_timegrid_mu_unst <- eta_timegrid_mu_unst -
        x$state$fitted_timegrid + fitted_timegrid
      eta_timegrid_mu <-
        (eta_timegrid_mu_unst - rep(long_bar, each = nsubj*n_w)) /
        rep(long_sds, each = nsubj*n_w)
      eta_timegrid_long <- rowSums(matrix(eta_timegrid_alpha *
                                            eta_timegrid_mu,
                                          nrow = nsubj*n_w,
                                          ncol = nmarker))
      eta_timegrid <- eta_timegrid_lambda + eta_timegrid_long
      eta_T_mu_unst <- eta_T_mu_unst - x$state$fitted_T + fitted_T
      eta_T_mu <- (eta_T_mu_unst - rep(long_bar, each = nsubj)) /
        rep(long_sds, each = nsubj)
      LogPost <- get_LogLik(eta_timegrid, eta_T_mu, eta) + x$prior(par)
      return(-1 * LogPost)
    }
    nu <- optimize(f = nu_update_opt, interval = c(0, 1))$minimum
  }
  b <- b + nu * Hs
  if ("random.effect" %in% class(x)) {
    b <- b - mean(b)
  }

  x$state$parameters[seq_len(b_p)] <- b
  x$state$fitted_timegrid <- drop(x$Xgrid %*% b)
  x$state$fitted.values <- drop(x$X %*% b)
  x$state$fitted_T <- drop(x$XT %*% b)
  x$state$hessian <- x_H
  x$state$edf <- sum_diag(x_H0 %*% Sigma)

  return(x$state)

}



# Updating sigma predictor ------------------------------------------------


update_mjm_sigma <- function(x, y, eta, eta_timegrid, eta_T_mu, survtime,
                             get_LogLik, update_nu, update_tau, edf,
                             iwls, ...) {

  b <- bamlss::get.state(x, "b")
  nu <- x$nu
  matrix_inv <- utils::getFromNamespace("matrix_inv", "bamlss")
  sum_diag <- utils::getFromNamespace("sum_diag", "bamlss")
  get.ic2 <- utils::getFromNamespace("get.ic2", "bamlss")
  tau2.optim <- utils::getFromNamespace("tau2.optim", "bamlss")

  if (iwls) {
    score <- drop(-1 + (y[[1]][, "obs"] - eta$mu)^2 / (exp(eta$sigma)^2))
    hess <- rep(2, nrow(y))
    z <- eta$sigma + 1/hess*score
  } else {
    x_score0 <- crossprod(x$X, -1 + (y[[1]][, "obs"] - eta$mu)^2 /
                            exp(eta$sigma)^2)
    x_H0 <- 2 * crossprod(x$X * drop((y[[1]][, "obs"] - eta$mu)/exp(eta$sigma)^2),
                          x$X * drop(y[[1]][, "obs"] - eta$mu))
  }

  # Updating Taus
  if(update_tau && (x$state$do.optim || !x$fixed)) {
    env <- new.env()
    par <- x$state$parameters
    edf0 <- edf - x$state$edf
    status <- attr(y, "status")

    # Function for Updating Taus
    # Defined here so that it accesses x_score and x_H
    tau_update_opt <- function(tau2) {
      par[x$pid$tau2] <- tau2
      if (iwls) {
        xhess <- crossprod(x$X *rep(2, nrow(y)), x$X)
        Sigma <- matrix_inv(1 * xhess, index = NULL)
      } else {
        x_score <- x_score0 + x$grad(score = NULL, par, full = FALSE)
        x_H <- x_H0 + x$hess(score = NULL, par, full = FALSE)
        Sigma <- matrix_inv(x_H, index = NULL)
        Hs <- Sigma %*% x_score
      }
      if(update_nu) {
        # Function for updating step length
        # Defined here so that it accesses Hs and par and fitted from parents
        nu_update_opt <- function(nu) {
          if (iwls) {
            b2 <- drop(Sigma %*% crossprod(x$X, z*2))
          } else {
            b2 <- drop(b + nu * Hs)
          }
          par[x$pid$b] <- b2
          fit <- drop(x$X %*% b2)
          eta$sigma <- eta$sigma - fitted(x$state) + fit
          LogPost <- get_LogLik(eta_timegrid, eta_T_mu, eta) + x$prior(par)
          return(-1 * LogPost)
        }
        nu <- optimize(f = nu_update_opt, interval = c(0, 1))$minimum
      }

      # Calculate the Information Criterion with given edf and nu
      if (iwls) {
        b2 <- drop(Sigma %*% crossprod(x$X, z*2))
        edf <- sum_diag((1 * xhess) %*% Sigma)
      } else {
        b2 <- drop(b + nu * Hs)
        edf1 <- sum_diag(x_H0 %*% Sigma)
      }
      fit <- drop(x$X %*% b2)
      eta$sigma <- eta$sigma - fitted(x$state) + fit
      edf <- edf0 + edf1
      logLik <- get_LogLik(eta_timegrid, eta_T_mu, eta)
      ic <- get.ic2(logLik, edf, n = length(eta$mu), type = "AICc")

      # First time running the function or new minimum
      if (is.null(env$ic_val) || (ic < env$ic_val)) {

        # Prepare the current state for output
        x$state$parameters[x$pid$b] <- b2
        x$state$parameters[x$pid$tau2] <- tau2
        x$state$fitted.values <- fit
        x$state$edf <- edf1
        x$state$hessian <- x_H

        assign("ic_val", ic, envir = env)
        assign("state", x$state, envir = env)
      }

      # For optimizer
      return(ic)
    }

    ic0 <- tau_update_opt(get.state(x, "tau2"))
    tau2 <- tau2.optim(tau_update_opt, start = get.state(x, "tau2"))
    return(env$state)

  }

  if (iwls) {
    xhess <- crossprod(x$X *rep(2, nrow(y)), x$X)
    Sigma <- matrix_inv(1 * xhess, index = NULL)
  } else {
    # Newton-Raphson
    x_score <- x_score0 + x$grad(score = NULL, x$state$parameters, full = FALSE)
    x_H <- x_H0 + x$hess(score = NULL, x$state$parameters, full = FALSE)
    Sigma <- matrix_inv(x_H, index = NULL)
    Hs <- Sigma %*% x_score
  }


  if(update_nu) {
    par <- x$state$parameters
    # Function for updating step length
    # Defined here so that it accesses Hs and par and fitted from parents
    nu_update_opt <- function(nu) {
      if(iwls) {
        b2 <- drop(Sigma %*% crossprod(x$X, z*2))
      } else {
        b2 <- drop(b + nu * Hs)
      }
      par[x$pid$b] <- b2
      fit <- drop(x$X %*% b2)
      eta$sigma <- eta$sigma - fitted(x$state) + fit
      LogPost <- get_LogLik(eta_timegrid, eta_T_mu, eta) + x$prior(par)
      return(-1 * LogPost)
    }
    nu <- optimize(f = nu_update_opt, interval = c(0, 1))$minimum
  }
  if (iwls) {
    b <- drop(Sigma %*% crossprod(x$X, z*2))
  } else {
    b <- b + nu * Hs
  }

  x$state$parameters[x$pid$b] <- b
  x$state$fitted.values <- drop(x$X %*% b)
  if (iwls) {
    x$state$hessian <- xhess
    x$state$edf <- sum_diag((1 * xhess) %*% Sigma)
  } else {
    x$state$hessian <- x_H
    x$state$edf <- sum_diag(x_H0 %*% Sigma)
  }


  return(x$state)

}
