
# MJM_mcmc ----------------------------------------------------------------

MJM_mcmc <- function(x, y, family, start = NULL, weights = NULL, offset = NULL,
                     n.iter = 1200, burnin = 200, thin = 1, step = 20,
                     nu_sampler = 1, verbose = FALSE, accthreshold = -1,
                     prop_list = NULL, std_surv = TRUE,
                     ...)
{

  # Set starting values for the sampling
  if(!is.null(start)) {
    if(is.matrix(start)) {
      if(any(i <- grepl("Mean", colnames(start))))
        start <- start[, i]
      else stop("the starting values should be a vector not a matrix!")
    }
    x <- bamlss::set.starting.values(x, start)
  }

  # Names of parameters/predictors and other attributes
  nx <- names(x)
  nmarker <- attr(y, "nmarker")
  take_last <- attr(y, "take_last")
  survtime <- y[[1]][, "time"][take_last]
  nsubj <- length(survtime)
  gq_weights <- attr(y, "gq_weights")
  nw <- length(gq_weights)
  status <- attr(y, "status")

  ## Compute additive predictors.
  eta <- bamlss:::get.eta(x, expand = FALSE)


  ## For the time dependent part, compute
  ## predictors based on the time grid.
  # mu
  eta_timegrid_mu <- 0
  eta_timegrid_mu_unst <- 0
  eta_T_mu <- 0
  eta_T_mu_unst <- 0
  if(length(x$mu$smooth.construct)) {
    for(j in names(x$mu$smooth.construct)) {
      b <- get.par(x$mu$smooth.construct[[j]]$state$parameters, "b")
      x$mu$smooth.construct[[j]]$state$fitted_timegrid <-
        drop(x$mu$smooth.construct[[j]]$Xgrid %*% b)
      x$mu$smooth.construct[[j]]$state$fitted_T <-
        drop(x$mu$smooth.construct[[j]]$XT %*% b)
      eta_timegrid_mu_unst <- eta_timegrid_mu_unst +
        x$mu$smooth.construct[[j]]$state$fitted_timegrid
      eta_T_mu_unst <- eta_T_mu_unst + x$mu$smooth.construct[[j]]$state$fitted_T
    }
  }
  # Calculate standardization constants
  if (std_surv) {
    x_std <- matrix(eta_timegrid_mu_unst, ncol = nmarker)
    long_bar <- colMeans(x_std)
    long_sds <- apply(x_std, 2, sd)
  } else {
    long_bar <- rep(0, nmarker)
    long_sds <- rep(1, nmarker)
  }
  attr(eta, "std_long") <- list(
    "long_bar" = long_bar,
    "long_sds" = long_sds
  )
  eta_timegrid_mu <-
    (eta_timegrid_mu_unst - rep(long_bar, each = nsubj * nw)) /
    rep(long_sds, each = nsubj * nw)
  eta_T_mu <- (eta_T_mu_unst - rep(long_bar, each = nsubj)) /
    rep(long_sds, each = nsubj)

  # alpha
  eta_timegrid_alpha <- 0
  if(length(x$alpha$smooth.construct)) {
    for(j in names(x$alpha$smooth.construct)) {
      b <- get.par(x$alpha$smooth.construct[[j]]$state$parameters, "b")
      x$alpha$smooth.construct[[j]]$state$fitted_timegrid <-
        drop(x$alpha$smooth.construct[[j]]$Xgrid %*% b)
      eta_timegrid_alpha <- eta_timegrid_alpha +
        x$alpha$smooth.construct[[j]]$state$fitted_timegrid
    }
  }
  # lambda
  eta_timegrid_lambda <- 0
  if(length(x$lambda$smooth.construct)) {
    for(j in names(x$lambda$smooth.construct)) {
      b <- get.par(x$lambda$smooth.construct[[j]]$state$parameters, "b")
      x$lambda$smooth.construct[[j]]$state$fitted_timegrid <-
        drop(x$lambda$smooth.construct[[j]]$Xgrid %*% b)
      eta_timegrid_lambda <- eta_timegrid_lambda +
        x$lambda$smooth.construct[[j]]$state$fitted_timegrid
    }
  }

  # Eta predictors
  eta_timegrid_long <- rowSums(matrix(eta_timegrid_alpha*eta_timegrid_mu,
                                      nrow = nsubj*nw, ncol = nmarker))
  eta_timegrid <- eta_timegrid_lambda + eta_timegrid_long
  eta_T_long <- rowSums(matrix(eta$alpha*eta_T_mu, nrow = nsubj,
                               ncol = nmarker))
  eta_T <- eta$lambda + eta$gamma + eta_T_long

  # Old logLikelihood and prior.
  sum_Lambda <- (survtime/2 * exp(eta$gamma)) %*%
    (diag(nsubj)%x%t(gq_weights))%*%
    exp(eta_timegrid)
  logLik_old <- drop(status %*% eta_T - sum_Lambda) +
    sum(dnorm(y[[1]][, "obs"], mean = eta$mu, sd = exp(eta$sigma),
              log = TRUE))

  # Fct for saving acceptance probability on the right scale
  transform_acceptprop <- function(x) {
    if(is.null(x)) return(0)
    if(is.na(x)) return(0)
    x <- exp(x)
    if(x < 0)
      x <- 0
    if(x > 1)
      x <- 1
    x
  }

  # Allow different steplenghts per model term
  for(j in names(x)) {
    for(sj in names(x[[j]]$smooth.construct)) {
      x[[j]]$smooth.construct[[sj]]$state$nu <- nu_sampler
      x[[j]]$smooth.construct[[sj]]$low_acc <- 0
    }
  }

  ## Process iterations
  if (burnin < 1) burnin <- 1
  if (burnin > n.iter) burnin <- floor(n.iter * 0.1)
  if (thin < 1) thin <- 1
  iterthin <- as.integer(seq(burnin, n.iter, by = thin))

  ## Samples.
  samps <- list()
  for (i in nx) {
    samps[[i]] <- list()
    for(j in names(x[[i]]$smooth.construct)) {
      samps[[i]][[j]] <- list(
        "samples" = matrix(NA, nrow = length(iterthin),
          ncol = length(x[[i]]$smooth.construct[[j]]$state$parameters)),
        "edf" = rep(NA, length = length(iterthin)),
        "alpha" = rep(NA, length = length(iterthin)),
        "accepted" = rep(NA, length = length(iterthin))
      )
      colnames(samps[[i]][[j]]$samples) <-
        names(x[[i]]$smooth.construct[[j]]$state$parameters)
    }
  }
  logLik.samps <- logPost.samps <- rep(NA, length = length(iterthin))

  nstep <- step
  step <- floor(n.iter / step)

  ptm <- proc.time()
  for(iter in 1:n.iter) {

    if (!is.null(prop_list)) {
      nx_iter <- 1
    }
    if(save <- iter %in% iterthin) {
      js <- which(iterthin == iter)
    }

    if(verbose) {
      cat("Iteration", iter, "\n")
    }

    for (i in nx) {

      if(!is.null(prop_list)) {
        j_iter <- 1
      }
      for (j in names(x[[i]]$smooth.construct)) {
        p_state <- propose_mjm(predictor = i,
                               x = x[[i]]$smooth.construct[[j]], y = y,
                               eta = eta, eta_timegrid = eta_timegrid,
                               eta_T = eta_T, eta_T_mu = eta_T_mu,
                               eta_T_mu_unst = eta_T_mu_unst,
                               eta_timegrid_alpha = eta_timegrid_alpha,
                               eta_timegrid_mu_unst = eta_timegrid_mu_unst,
                               eta_timegrid_mu = eta_timegrid_mu,
                               eta_timegrid_long = eta_timegrid_long,
                               eta_timegrid_lambda = eta_timegrid_lambda,
                               survtime = survtime, logLik_old = logLik_old,
                               nsubj = nsubj, gq_weights = gq_weights,
                               status = status, nmarker = nmarker,
                               verbose_sampler = verbose,
                               prop = if(!is.null(prop_list)) {
                                 prop_list[[iter]][[nx_iter]][[j_iter]]
                               } else NULL)

        # If accepted, set current state to proposed state
        accepted <- if(!is.null(prop_list)) {
          attr(prop_list[[iter]][[nx_iter]][[j_iter]], "acc")
        } else {if(is.na(p_state$xstate$alpha)){
          FALSE
        } else {
          log(runif(1)) <= p_state$xstate$alpha
        }}

        # Check whether the step-size should be decreased
        if(is.na(p_state$xstate$alpha) ||
           exp(p_state$xstate$alpha) < accthreshold) {
          x[[i]]$smooth.construct[[j]]$low_acc <-
            x[[i]]$smooth.construct[[j]]$low_acc + 1
        } else {
          x[[i]]$smooth.construct[[j]]$low_acc <- 0
        }
        if (x[[i]]$smooth.construct[[j]]$low_acc == 10) {
          x[[i]]$smooth.construct[[j]]$state$nu <-
            x[[i]]$smooth.construct[[j]]$state$nu / 10
          x[[i]]$smooth.construct[[j]]$low_acc <- 0
          cat("HAPPENED\n")
        }

        if (accepted) {

          # Update the etas
          switch(i, "lambda" = {
            eta_T <- p_state$etas$eta_T
            eta_timegrid <- p_state$etas$eta_timegrid
            eta_timegrid_lambda <- p_state$etas$eta_timegrid_lambda
          }, "gamma" = {
            eta_T <- p_state$etas$eta_T
          }, "alpha" = {
            eta_T <- p_state$etas$eta_T
            eta_timegrid <- p_state$etas$eta_timegrid
            eta_timegrid_long <- p_state$etas$eta_timegrid_long
            eta_timegrid_alpha <- p_state$etas$eta_timegrid_alpha
          }, "mu" = {
            eta_T <- p_state$etas$eta_T
            eta_T_mu <- p_state$etas$eta_T_mu
            eta_T_mu_unst <- p_state$etas$eta_T_mu_unst
            eta_timegrid <- p_state$etas$eta_timegrid
            eta_timegrid_long <- p_state$etas$eta_timegrid_long
            eta_timegrid_mu <- p_state$etas$eta_timegrid_mu
            eta_timegrid_mu_unst <- p_state$etas$eta_timegrid_mu_unst
          })
          eta <- p_state$etas$eta

          # Update likelihood and state
          logLik_old <- p_state$logLik
          x[[i]]$smooth.construct[[j]]$state <- p_state$xstate
        }


        ## Save the samples and acceptance.
        if(save) {
          samps[[i]][[j]]$samples[js, ] <-
            x[[i]]$smooth.construct[[j]]$state$parameters
          samps[[i]][[j]]$edf[js] <- x[[i]]$smooth.construct[[j]]$state$edf
          samps[[i]][[j]]$alpha[js] <-
            transform_acceptprop(p_state$xstate$alpha)
          samps[[i]][[j]]$accepted[js] <- accepted
        }
        if (!is.null(prop_list)){
          j_iter <- j_iter + 1
        }
      }

      if(!is.null(prop_list)) {
        nx_iter <- nx_iter + 1
      }
    }
    if(save) {
      logLik.samps[js] <- logLik_old
      logPost.samps[js] <- as.numeric(logLik.samps[js] +
                                        bamlss:::get.log.prior(x))
    }
  }

  for(i in names(samps)) {
    for(j in names(samps[[i]])) {
      samps[[i]][[j]] <- do.call("cbind", samps[[i]][[j]])
      cn <- if(j == "model.matrix") {
        paste(i, "p", j, colnames(samps[[i]][[j]]), sep = ".")
      } else {
        paste(i, "s", j, colnames(samps[[i]][[j]]), sep = ".")
      }
      colnames(samps[[i]][[j]]) <- cn
    }
    samps[[i]] <- do.call("cbind", samps[[i]])
  }
  samps$logLik <- logLik.samps
  samps$logPost <- logPost.samps
  samps <- do.call("cbind", samps)

  ## Compute DIC. #
  dev <- -2 * logLik.samps
  mpar <- apply(samps, 2, mean, na.rm = TRUE)
  names(mpar) <- colnames(samps)
  ll <- family$p2logLik(mpar)
  mdev <- -2 * ll
  pd <- mean(dev) - mdev
  DIC <- mdev + 2 * pd
  samps <- cbind(samps,
                 "DIC" = rep(DIC, length.out = nrow(samps)),
                 "pd" = rep(pd, length.out = nrow(samps))
  )
  if (std_surv) {

    # Scaling of the survival design matrices
    # NOTE: The samples are adapted so that the parameters are on the original
    # scale. This makes interpretation and prediction easier. Parameters from
    # the optimizing step, however, are not yet adjusted as they are needed as
    # starting values for the MCMC algorithm. Note that the implementation so
    # far focusses only on one gamma covariate without smooth constructor term.
    # If a more complex model is fit with JMbamlss, please adjust the following
    # code accordingly.
    # In general: Scale the parameters beta by the standard deviation sd.
    # Subtract beta*bar/sd from the intercept.

    # Rescale the alpha samples (accounting for standardization)
    alpha_cols <- grep("alpha.+marker", colnames(samps))
    for (i in seq_along(alpha_cols)) {
      samps[, alpha_cols[i]] <- samps[, alpha_cols[i]] / long_sds[i]
    }

    # Calculate gamma intercept adaption for standardization
    surv_int <- apply(samps[, alpha_cols], 1, "%*%", long_bar)

    # Rescale the gamma covariates samples
    gamma_bar <- x$gamma$smooth.construct$model.matrix$w_bar
    gamma_sd <- x$gamma$smooth.construct$model.matrix$w_sd
    if (length(gamma_bar) > 1)
      warning("Standardization is not tested for this case of gamma.")
    gamma_names <- colnames(x$gamma$smooth.construct$model.matrix$X)[-1]
    for (j in seq_along(gamma_names)) {
      gamma_col <- grep(paste0("gamma.+", gamma_names[j]), colnames(samps))
      samps[, gamma_col] <- samps[, gamma_col] / gamma_sd[j]
      surv_int <- surv_int + samps[, gamma_col] * gamma_bar[j]
    }

    # Adapt the intercept
    col_gamma_int <- grep("gamma\\.p\\.model\\.matrix\\.\\(Intercept",
                          colnames(samps))
    samps[, col_gamma_int] <- samps[, col_gamma_int] - surv_int

  }
  samps[is.na(samps)] <- 0


  return(as.mcmc(samps))
}

