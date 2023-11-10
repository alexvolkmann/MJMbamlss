

#' New Simulation Function For Multivariate JMs Based On FPCs
#'
#' Adapt the structure given by simJM function in bamlss.
#'
#' @param nsub Number of subjects.
#' @param times Vector of time points.
#' @param probmiss Probability of missingness.
#' @param max_obs Maximal number of observations per individual and marker.
#'   Defaults to no upper limit.
#' @param maxfac Factor changing the uniform censoring interval.
#' @param nmark Number of markers.
#' @param long_assoc Longitudinal association between the markers (Defaults to
#'   "FPC"). If "splines" or "param", then specify the normal covariance matrix
#'   with argument 're_cov_mat' and include the random effects in argument mu.
#'    If "FPC", then principal components are used to model the association
#'    structure.
#' @param M Number of principal components.
#' @param FPC_bases FunData object. If supplied, use the contained FPC as basis
#'   for the association structure.
#' @param FPC_evals Vector of eigenvalues. If supplied, use the provided
#'   eigenvalues for the association structure.
#' @param mfpc_args List containing the named arguments "type", "eFunType",
#'   "ignoreDeg", "eValType" of function simMultiFunData and "eValScale" for
#'   scaling the eigenvalues.
#' @param re_cov_mat If supplied, a covariance matrix to use for drawing the
#'   random effects needed for the association structure.
#' @param ncovar Number of covariates.
#' @param lambda Additive predictor of time-varying survival covariates.
#' @param gamma Additive predictor of time-constant survival covariates.
#' @param alpha List of length nmark containing the additive predictors of the
#'   association.
#' @param mu List of length nmark containing the additive predictors of the
#'   longitudinal part.
#' @param sigma Additive predictor of the variance.
#' @param tmax Maximal time point of observations.
#' @param seed Seed for reproducibility.
#' @param mfpc_args List containing the named arguments "type", "eFunType",
#'   "ignoreDeg", "eValType" of function simMultiFunData and "eValScale" for
#'   scaling the eigenvalues.
#' @param full Create a wide-format data.frame and a short one containing only
#'   survival info.
#' @param file Name of the data file the generated data set should be stored
#'   into (e.g., "simdata.RData") or NULL if the dataset should directly be
#'   returned in R.
#' @export
simMultiJM <- function(nsub = 300, times = seq(0, 120, 1), probmiss = 0.75,
                       max_obs = length(times), maxfac = 1.5, nmark = 2,
                       long_assoc = c("FPC", "splines", "param"), M = 6,
                       FPC_bases = NULL, FPC_evals = NULL,
                       mfpc_args = list(type = "split", eFunType = "Poly",
                                        ignoreDeg = NULL, eValType = "linear",
                                        eValScale = 1),
                       re_cov_mat = NULL, ncovar = 2,
                       lambda = function(t, x) {
                         1.4*log((t + 10)/1000)
                       },
                       gamma = function(x) {
                         - 1.5 + 0.3*x[, 1]
                       },
                       alpha = rep(list(function(t, x) {
                         0.3 + 0*t
                       }), nmark),
                       mu = rep(list(function(t, x){
                         1.25 + 0.6*sin(x[, 2]) + (-0.01)*t
                       }), nmark),
                       sigma = function(t, x) {
                         0.3 + 0*t + I(x$marker == "m2")*0.2
                       },
                       tmax = NULL, seed = NULL,
                       full = FALSE, file = NULL){

  long_assoc <- match.arg(long_assoc)
  # Some basic checks
  if(length(alpha) != length(mu)) {
    stop("alpha and mu must have same length.\n")
  }
  if(length(mu) != nmark) {
    stop("Predictors must be specified for all markers.\n")
  }
  if(!is.null(re_cov_mat) & long_assoc == "FPC") {
    stop("Either REs or PCREs need to be specified properly.\n")
  }
  if(long_assoc != "FPC" & is.null(re_cov_mat)) {
    stop("Specify the RE covariance matrix.\n")
  }
  if(is.null(tmax)){
    tmax <- max(times)
  }
  if(!is.null(FPC_bases)) {
    if (M < funData::nObs(FPC_bases)) {
      fpc_base <- funData::extractObs(FPC_bases, seq_len(M))
    } else {
      M <- funData::nObs(FPC_bases)
    }
    if (!is.null(FPC_evals)) {
      if (length(FPC_evals) < funData::nObs(FPC_bases)) {
        stop("Wrong number of eigenvalues supplied.")
      } else {
        FPC_evals <- FPC_evals[seq_len(M)]
      }
    }
  }


  ## specify censoring function (aus CoxFlexBoost) changed into uniformly (as
  ## too much censoring)
  ## added censoring at tmax
  censoring <- function(time, tmax, maxfac){
    ## censoring times are independent uniformly distributed
    censor_time <- runif(n = length(time), min = 0, max = maxfac*tmax)
    censor_time <- ifelse(censor_time > tmax, tmax, censor_time)
    event <- (time <= censor_time)
    survtime <- apply(cbind(time, censor_time), 1, min)
    ## return matrix of observed survival times and event indicator
    return(cbind(survtime, event))
  }


  ## introduce random missings
  ## (excluding the first measurement to ensure all subjects remain in the data
  # but if there is a maximal number of observations than we don't necessarily
  # have the first longitudinal observation)
  miss_fct <- function(data, prop, max_obs, obstime = "obstime"){

    select <- which(data[[obstime]] > 0)
    n <- length(select)
    n_miss <- round(prop*n, 0)
    miss <- sample(select, n_miss)
    if (length(miss) > 0) {
      data <- data[-miss,]
    }
    if (max_obs < length(times)) {
      id_mark <- split(data, interaction(data$id, data$marker))
      select <- do.call(c, lapply(id_mark, function(x) {
        if (nrow(x) < max_obs) {
          rep(TRUE, times = nrow(x))
        } else {
          c(TRUE, sample(c(rep(TRUE, times = max_obs - 1),
                           rep(FALSE, times = nrow(x) - max_obs)),
                         nrow(x) - 1, replace = FALSE))
        }
      }))
      data <- data[select, ]
    }
    return(data)
  }


  ## generate baseline covariates
  gen_x <- function(nsub, ncovar){
    x <- matrix(data = NA, nrow = nsub, ncol = ncovar)
    for(i in seq_len(ncovar)) {
      x[, i] <- runif(n = nsub, min = -3, max = 3)
    }
    colnames(x) <- paste0("x", seq_len(ncovar))
    facs <- list()
    for(i in seq_len(ncovar)) {
      facs[[i]] <- sample(c(0, 1), size = nsub, replace = TRUE)
    }
    names(facs) <- paste0("x", ncovar + seq_len(ncovar))
    data.frame(x, facs)
  }

  ## generate parametric random effects
  gen_b <- function (nsub, re_cov_mat) {

    r <- mvtnorm::rmvnorm(n = nsub, sigma = re_cov_mat)
    colnames(r) <- paste0("r", seq_len(nrow(re_cov_mat)))
    r <- data.frame(r)
    out <- list(r, list())
  }



  ## generate the multivariate functional principal component basis
  mfpc <- function(argvals, mfpc_args, M) {
    bases <- switch(mfpc_args$type,
      split = simMultiSplit(argvals = argvals, M = M,
                            eFunType = mfpc_args$eFunType,
                            ignoreDeg = mfpc_args$ignoreDeg,
                            eValType = mfpc_args$eValType,
                            s = mfpc_args$mfpc_seed),
      weighted = simMultiWeight(argvals = argvals, M = M,
                                eFunType = mfpc_args$eFunType,
                                ignoreDeg = mfpc_args$ignoreDeg,
                                eValType = mfpc_args$eValType,
                                alpha = mfpc_args$mfpc_seed),
      stop(paste("Choose either 'split' or 'weighted' for the simulation",
                 "of multivariate functional data.")))
    if (M == 1) {
      bases <- funData::multiFunData(lapply(bases, "/",
                                            sqrt(funData::norm(bases))))
    }
    bases
  }

  ## generate functional principal component based random effects
  gen_fpc <- function(times, nsub, M, FPC_bases = NULL, FPC_evals = NULL,
                      mfpc_args, tmax, seed = NULL){
    if(!is.null(seed)) set.seed(seed)

    if (!is.null(FPC_evals)) {
      evals <- FPC_evals
    } else {
      evals <- funData::eVal(M = M, type = mfpc_args$eValType)
      evals <- mfpc_args$eValScale * evals
    }

    scores <- mvtnorm::rmvnorm(nsub, sigma = diag(evals),
                               method="chol")
    colnames(scores) <- paste0("s", 1:M)

    if (!is.null(FPC_bases)) {
      b_set <- list("FPC_bases" = FPC_bases, "M" = M, "evals" = evals)
    } else {
      mfpc_seed <- switch(mfpc_args$type,
                          "split" = sample(c(-1, 1), nmark, 0.5),
                          "weight" = stats::runif(nmark, 0.2, 0.8))
      b_set <- c(list(tmin = min(c(times, tmax)), tmax = tmax, M = M,
                      mfpc_seed = mfpc_seed, evals = evals), mfpc_args)
    }
    return(list(scores, b_set))
  }

  ## compute predictors
  ## individual longitudinal trajectories
  mu_fun <-  function(time, x, r, mu, b_set){

    # parametric random effects
    if (length(b_set) == 0){
      out <- lapply(mu, function (mu_k) {
        mu_k(t = time, x = x, r = r)
      })

    } else if (long_assoc == "FPC") {

      # duplicate scores for the multiple integration points
      if(is.null(dim(r))){
        r <- matrix(r, nrow = length(time), ncol = b_set$M, byrow=TRUE)
      }

      # FPCs are normalized on the interval up to tmax so do a 'last value
      # carried forward' type of evaluation
      time_out <- ifelse(time > tmax, tmax, time)

      # Evaluate the functional principal component bases for different markers
      if(is.null(b_set$FPC_bases)) {
        pc_bases <- lapply(
          mfpc(argvals = rep(list(c(b_set$tmin, time_out, b_set$tmax)),
                             length(mu)),
               mfpc_args = b_set, M = b_set$M),
          function (fundat){
            t(fundat@X)[-c(1, 2+length(time)), ]
          })
      } else {
        pc_bases <- lapply(b_set$FPC_bases, eval_fundata,
                           evalpoints = time_out)
      }

      out <- mapply(function (mu_k, pc_k) {
        mu_k(t = time, x = x) + apply(pc_k*r, 1, sum)
      }, mu_k = mu, pc_k = pc_bases, SIMPLIFY = FALSE)

    } else {

      # Use last-value-carried forward if time is larger than tmax otherwise
      # error in splineDesign
      time_out <- ifelse(time > tmax, tmax, time)

      # Slightly different for spline based random effects
      spline_base <- splines::splineDesign(
        knots = b_set$knots,
        x = time_out, ord = 4, derivs = 0)

      # random effects are different over the markers but the splines are the
      # same
      # duplicate scores for the multiple integration points
      if(dim(r)[1] == 1){
        r <- matrix(as.matrix(r), nrow = length(time), ncol = ncol(r), byrow=TRUE)
      }
      random_list <- lapply(seq_len(nmark), function(ind) {
        r[, (ind-1)*ncol(spline_base) + seq_len(ncol(spline_base))]
      })

      out <- mapply(function (mu_k, ran_k) {
        mu_k(t = time, x = x) + apply(spline_base * ran_k, 1, sum)
      }, mu_k = mu, ran_k = random_list, SIMPLIFY = FALSE)
    }

    out

  }

  ## full hazard
  hazard <-  function(time, x, r, ...){

    mu_list <- mu_fun(time, x, r, mu, b_set)

    exp(lambda(time, x) + gamma(x) +
          Reduce('+', mapply(function(alpha_k, mu_k) {
            alpha_k(time, x)*mu_k
          }, alpha_k = alpha, mu_k = mu_list,
          SIMPLIFY = FALSE)))
  }


  ## Code from bamlss package, slightly adapted
  rJM <- function (hazard, censoring, x, r, subdivisions = 1000, tmin = 0,
            tmax, maxfac, file = NULL, ...)
  {
    nsub <- nrow(x)
    time <- rep(NA, nsub)
    Hazard <- function(hazard, time, x, r) {
      integrate(hazard, 0, time, x = x, r = r,
                subdivisions = subdivisions, stop.on.error = FALSE)$value
    }
    InvHazard <- function(Hazard, hazard, x, r, tmin, tmax) {
      negLogU <- -log(runif(1, 0, 1))
      rootfct <- function(time) {
        negLogU - Hazard(hazard, time, x, r)
      }
      if (rootfct(tmin) < 0) {
        return(0)
      }
      else {
        root <- try(uniroot(rootfct, interval = c(0, tmax))$root,
                    silent = TRUE)
        root <- if (inherits(root, "try-error")) {
          tmax + 0.01
        }
        else {
          root
        }
      }
      return(root)
    }
    cumhaz <- rep(NA, nsub)
    survprob <- rep(NA, nsub)
    for (i in 1:nsub) {
      time[i] <- InvHazard(Hazard, hazard, x[i, ], r[i, ],
                           tmin, tmax)
      cumhaz[i] <- Hazard(hazard, time[i], x[i, ], r[i, ])
      survprob[i] <- exp((-1) * cumhaz[i])
    }
    time_event <- censoring(time, tmax, maxfac)
    haz <- rep(NA, nsub)
    for (i in 1:nsub){
      haz[i] <- hazard(time = time_event[i, 1], x[i, ], r[i, ])
    }
    data_short <- data.frame(survtime = time_event[, 1],
                             event = time_event[, 2], x, r, cumhaz = cumhaz,
                             hazard = haz)
    names(data_short) <- gsub(".", "", names(data_short), fixed = TRUE)
    return(data_short)
  }


  # generate input
  id <- rep(1:nsub, each=length(times))
  if(!is.null(seed)){
    set.seed(seed)
  }
  x <- gen_x(nsub, ncovar = ncovar)


  if (long_assoc != "FPC") {
    temp <- gen_b(nsub = nsub, re_cov_mat = re_cov_mat)
    if (long_assoc == "splines") {

      nk <- nrow(re_cov_mat) / nmark - 2
      x_r <- tmax - times[1]
      xl <- times[1] - x_r*0.001
      xu <- tmax + x_r*0.001
      dx <- (xu - xl)/(nk - 1)
      kn <- seq(xl - dx * 3, xu + dx * 3, length = nk + 2 * 3)

      temp[[2]] <- list(knots = kn,
                        tmin = times[1])
    }
  } else {
    temp <- gen_fpc(times = times, nsub = nsub, M = M, FPC_bases = FPC_bases,
                    FPC_evals = FPC_evals, mfpc_args = mfpc_args, tmax = tmax)
  }
  r <- temp[[1]]
  b_set <- temp[[2]]

  data_short <- rJM(hazard, censoring, x, r, tmin = times[1], tmax = tmax,
                    maxfac = maxfac)


  ## Create the full simulated data
  data_base <- cbind(id, data_short[id,], obstime = rep(times, nsub))
  i <- !duplicated(data_base$id)
  data_base$id <- as.factor(data_base$id)

  # gamma and lambda have only joint intercept which is estimated in
  # predictor gamma
  f_lambda <- lambda(data_short$survtime)
  data_base$lambda <- lambda(data_base$obstime) - mean(f_lambda)
  data_base$gamma <- gamma(x[id, ]) + mean(f_lambda)

  data_long <- do.call(rbind, rep(list(data_base), nmark))
  data_long$marker <- factor(rep(paste0("m", seq_len(nmark)),
                                 each = length(id)))
  data_long$mu <- do.call(c, mu_fun(data_base$obstime, x[id, ], r[id,], mu,
                                    b_set))
  data_long$alpha <- do.call(c, lapply(alpha, function(alpha_k) {
                       alpha_k(data_base$obstime, x[id, ])
                     }))
  data_long$sigma <- sigma(t = data_long$obstime, x = data_long)
  if (long_assoc == "FPC") {
    if (is.null(FPC_bases)) {
      fpcs <- do.call(rbind, lapply(mfpc(argvals = rep(list(data_base$obstime),
                                                       nmark),
                                         mfpc_args = b_set, M = M),
                                    function (mark) {
                                      t(mark@X)
                                    }))
    } else {
      fpcs <- do.call(rbind, lapply(b_set$FPC_bases, eval_fundata,
                                    evalpoints = data_base$obstime))
    }

    data_long <- cbind(data_long,
                       fpc = fpcs,
                       wfpc = t(t(fpcs)*sqrt(b_set$evals)))
  } else if (long_assoc == "splines") {

    # add spline based random effects
    spline_base <- splines::splineDesign(
      knots = b_set$knots,
      x = data_base$obstime, ord = 4, derivs = 0)
    random_list <- lapply(seq_len(nmark), function(ind) {
      r[id, (ind-1)*ncol(spline_base) + seq_len(ncol(spline_base))]
    })

    data_long$fre <- do.call(c, lapply(random_list, function (ran_k) {
      apply(spline_base * ran_k, 1, sum)
    }))

  }

  # Hypothetical longitudinal data
  data_hypo <- data_long

  # Data at survival time
  data_short$id <- factor(seq_len(nsub))
  data_short$lambda <- lambda(data_short$survtime) - mean(f_lambda)
  data_short$gamma <- gamma(x) + mean(f_lambda)
  shortmu <- do.call(c, mu_fun(data_short$survtime, x, r, mu, b_set))
  shortalpha <- do.call(c, lapply(alpha, function(alpha_k) {
    alpha_k(data_short$survtime, x)
  }))
  data_short <- do.call(rbind, rep(list(data_short), nmark))
  data_short$marker <- factor(rep(paste0("m", seq_len(nmark)),
                                  each = nrow(x)))
  data_short$alpha <- shortalpha
  data_short$mu <- shortmu
  data_short$sigma <- sigma(t = data_short$survtime, x = data_short)
  attr(data_short, "f_lambda") <- f_lambda

  # censoring
  data_long <- data_long[data_long$obstime <= data_long$survtime,]

  # saving data without longitudinal missings
  data_full <- data_long

  # inducing longitudinal missings
  data_long <- miss_fct(data_long, probmiss, max_obs)

  # Draw longitudinal observations
  data_long$y <- rnorm(nrow(data_long), data_long$mu, sd = exp(data_long$sigma))


  if(full){

    if (long_assoc == "FPC") {
      # FPC basis functions
      if (is.null(FPC_bases)) {
        FPC_bases <- mfpc(argvals = rep(list(times), nmark), mfpc_args = b_set,
                          M = M)
      }
    }

    d <- list(data = data_long, data_full = data_full, data_hypo = data_hypo,
              fpc_base = if (long_assoc == "FPC") FPC_bases,
              data_short = data_short)
  } else {
    d <- data_long
  }
  if(!is.null(file)) {
    save(d, file = file)
    invisible(d)
  } else {
    return(d)
  }
}




# Function to create a MFPCA basis using splits ---------------------------


## Code from Clara Happ-Kurz' package funData and slightly adapted
simMultiSplit <- function (argvals, M, eFunType, ignoreDeg = NULL, eValType,
                           s) {
  if (any(c(length(M), length(eFunType), length(eValType)) != 1))
    stop("argvals, M, eFunType, eValType must all be of length 1!")
  p <- length(argvals)
  x <- vector("list", length = length(argvals))
  splitVals <- rep(NA, length(argvals) + 1)
  x[[1]] <- unlist(argvals[[1]])
  splitVals[1:2] <- c(0, length(x[[1]]))
  if (p > 1) {
    for (i in 2:p) {
      x[[i]] <- unlist(argvals[[i]])
      x[[i]] <- argvals[[i]] - min(argvals[[i]]) + max(x[[i - 1]])
      splitVals[i + 1] <- splitVals[i] + length(x[[i]])
    }
  }
  f <- funData::eFun(unlist(x), M, ignoreDeg = ignoreDeg, type = eFunType)
  trueFuns <- vector("list", p)
  for (j in seq_len(p)) trueFuns[[j]] <- funData::funData(
    argvals[[j]], s[j] * f@X[, (1 + splitVals[j]):splitVals[j + 1],
                             drop = FALSE])
  return(funData::multiFunData(trueFuns))
}


# Function to crate a MFPCA basis using weighting -------------------------


## Code from Clara Happ-Kurz' package funData and slightly adapted
simMultiWeight <- function (argvals, M, eFunType, ignoreDeg = NULL, eValType,
                            alpha) {
  p <- length(argvals)
  dimsSupp <- foreach::foreach(j = seq_len(p), .combine = "c") %do% {
    length(argvals[[j]])
  }
  if (any(dimsSupp > 2)) {
    stop(paste0("Function simMultiWeight: method is not implemented for ",
                "objects of dimension > 2!"))
  }
  if (p > 1) {
    if (isTRUE(do.call(all.equal, lapply(M, prod)))) {
      Mtotal <- prod(M[[1]])
    }
    else stop("Function simMultiWeight: basis dimensions must be equal!")
  }
  else {
    Mtotal <- prod(M[[1]])
  }

  weight <- sqrt(alpha/sum(alpha))
  basis <- vector("list", p)
  for (j in seq_len(p)) {
    if (dimsSupp[j] == 1)
      basis[[j]] <- weight[j] * funData::eFun(argvals[[j]][[1]], M = M[[j]],
                                              ignoreDeg = ignoreDeg[[j]],
                                              type = eFunType[[j]])
    else basis[[j]] <- weight[j] * funData::tensorProduct(
      funData::eFun(argvals[[j]][[1]], M = M[[j]][1],
                    ignoreDeg = ignoreDeg[[j]][[1]], type = eFunType[[j]][1]),
      funData::eFun(argvals[[j]][[2]], M = M[[j]][2],
                    ignoreDeg = ignoreDeg[[j]][[2]],
                    type = eFunType[[j]][2]))
  }
  return(funData::multiFunData(basis))
}
