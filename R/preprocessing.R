
# Preprocessing Steps for Model Fitting -----------------------------------



# MFPCA object from data --------------------------------------------------

#' Preprocessing step to create MFPCA object
#'
#' This function takes the data und uses the residuals of marker-specific
#' additive models to estimate the covariance structure for a MFPCA
#'
#' @param data Data.frame such as returned by function simMultiJM.
#' @param uni_mean String to crate a formula for the univariate addtive models.
#' @param time String giving the name of the longitudinal time variable.
#' @param id String giving the name of the identifier.
#' @param marker String giving the name of the longitudinal marker variable.
#' @param M Number of mFPCs to compute in the MFPCA. If not supplied, it
#'  defaults to the maximum number of computable mFPCs.
#' @param weights TRUE if inverse sum of univariate eigenvals should be used as
#'  weights in the scalar product of the MFPCA. Defaults to FALSE (weights 1).
#' @param remove_obs Minimal number of observations per individual and marker to
#'  be included in the FPC estimation. Defaults to NULL (all observations). Not
#'  removing observations can lead to problems if the univariate variance
#'  estimate is negative and has to be truncated, then the scores for IDs with
#'  few observations cannot be estimated.
#' @param method Which package to use for the univariate FPCA. Either function
#'  adapted function 'fpca', 'FPCA' from package \code{fdapace}, 'fpca.sc' from
#'  package \code{refund}, or function 'PACE' from package \code{MFPCA}.
#' @param nbasis Number of B-spline basis functions for mean estimate for
#'  methods fpca, fpca.sc, PACE. For fpca.sc, PACE also bivariate smoothing of
#'  covariance estimate.
#' @param nbasis_cov Number of basis functions used for the bivariate
#'  smoothing of the covariance surface for method fpca.
#' @param bs_cov Type of spline for the bivariate smoothing of the covariance
#'  surface for method fpca. Default is symmetric fast covariance smoothing
#'  proposed by Cederbaum.
#' @param npc Number of univariate principal components to use in fpca.sc, PACE.
#' @param fve_uni Fraction of univariate variance explained for method FPCA.
#' @param pve_uni Proportion of univariate variance explained for methods
#'  fpca, fpca.sc, PACE.
#' @param fit MFPCA argument to return a truncated KL fit to the data. Defaults
#'  to FALSE.
#' @param max_time If supplied, forces the evaluation of the MFPCs up to maxtime.
#'  Only implemented for method = 'fpca'.
#' @param save_uniFPCA TRUE to attach list of univariate FPCAs as attribute to
#'  output. Defaults to FALSE.
#' @param save_uniGAM TRUE to attach list of univariate additive models used to
#'  calculate the residuals. Defaults to FALSE.
#' @export
#' @returns An object of class \code{MFPCAfit} with additional attributes
#'  depending on the arguments \code{save_uniFPCA}, \code{save_uniGAM},
#'  \code{fit}.
#' @examples
#' data(pbc_subset)
#' mfpca <- preproc_MFPCA(pbc_subset, uni_mean = paste0(
#'     "logy ~ 1 + sex + drug + s(obstime, k = 10, bs = 'ps') + ",
#'     "s(age, k = 10, bs = 'ps')"),
#'     pve_uni = 0.99, nbasis = 5, weights = TRUE, save_uniFPCA = TRUE)
preproc_MFPCA <- function (data, uni_mean = "y ~ s(obstime) + s(x2)",
                           time = "obstime", id = "id", marker = "marker",
                           M = NULL, weights = FALSE, remove_obs = NULL,
                           method = c("fpca", "fpca.sc", "FPCA", "PACE"),
                           nbasis = 10, nbasis_cov = nbasis, bs_cov = "symm",
                           npc = NULL, fve_uni = 0.99, pve_uni = 0.99,
                           fit = FALSE, max_time = NULL,
                           save_uniFPCA = FALSE, save_uniGAM = FALSE) {

  method <- match.arg(method)

  if (!is.null(remove_obs)) {
    few_obs <- apply(table(data[, id], data[, marker]), 1,
                     function (x) any(x < remove_obs))
    data <- droplevels(data[data[, id] %in% paste(which(!few_obs)), ])

  }

  marker_dat <- split(data, data$marker)


  # Check whether enough observations are available on each marker to be able
  # to estimate the same full interval on all markers
  maxtime <- sapply(marker_dat, function(x) max(x[, time]))
  if (!is.null(max_time)) {
    maxtime <- c(maxtime, max_time)
  }

  if (length(unique(maxtime)) > 1) {
    if (method != "fpca") {
      stop("Estimation of MFPCA for different univariate time intervals ",
           "has only been tested for method 'fpca'. Note that the design ",
           "matrices constructed in the joint model will not be correct.")
    } else {
      warning("Estimation of univariate FPCAs has to extrapolate for ",
              paste(paste0("marker", which(maxtime != max(maxtime))),
                    sep = ", "),
              " from ",
              paste(maxtime[which(maxtime != max(maxtime))], sep = ","),
              " to ", max(maxtime), ". Please check if appropriate.")
    }
  }
  argvals_pred <- lapply(marker_dat, function (x) {
    sort(unique(c(x[, time], max(maxtime))))
  })

  uni_mean <- as.formula(uni_mean)

  marker_mod <- lapply(marker_dat, function (mark) {
    bam(formula = uni_mean, data = mark)
  })
  marker_dat <- mapply(function (mark, mod) {
    mark$res <- mod$residuals
    mark
  }, mark = marker_dat, mod = marker_mod, SIMPLIFY = FALSE)

  if (method == "fpca.sc" | method == "fpca") {

    # Construct objects for fpca.sc function
    lY <- lapply(marker_dat, function (mark) {
      data.frame(".id" = mark[, id],
                 ".index" = mark[, time],
                 ".value" = mark$res)
    })

    # FPCA for each marker
    if (method == "fpca.sc") {
      FPCA <- lapply(lY, function(y) {
        refund::fpca.sc(ydata = y, pve = pve_uni, nbasis = nbasis, npc = npc,
                        var = TRUE)
      })
    } else {
      # For loop as it is more memory efficient
      FPCA <- list()
      for (m in seq_along(marker_dat)) {
        FPCA[[m]] <- fpca(ydata = lY[[m]], pve = pve_uni, nbasis = nbasis,
                          nbasis_cov = nbasis_cov,
                          bs_cov = bs_cov, npc = npc,
                          argvals_pred = argvals_pred[[m]])
      }
    }

    # Construct multivariate FunData and estimated FPCs
    mFData <- funData::multiFunData(lapply(FPCA, function (mark) {
      funData::funData(argvals = mark$argvals, X = mark$Yhat)
    }))
    uniExpansions <- lapply(FPCA, function (mark) {
      list(type = "given", functions = funData::funData(
        argvals = mark$argvals, X = t(mark$efunctions)))
    })
    if (is.null(M)) {
      M <- sum(sapply(FPCA, "[[", "npc"))
    }

    # Extract weights
    weight_vec <- 1/sapply(lapply(FPCA, "[[", "evalues"), sum)


  } else if (method == "FPCA") {

    # Construct objects for FPCA function
    ly <- lapply(marker_dat, function (mark) {
      mark <- mark[order(mark[, time]), ]
      split(mark$res, mark[, id])
    })
    lt <- lapply(marker_dat, function (mark) {
      mark <- mark[order(mark[, time]), ]
      split(mark[, time], mark[, id])
    })

    # FPCA for each marker
    FPCA <- mapply(fdapace::FPCA, Ly = ly, Lt = lt, SIMPLIFY = FALSE,
                   MoreArgs = list(optns = list(FVEthreshold = fve_uni,
                                                nRegGrid = 101)))

    # Construct multivariate FunData and estimated FPCs
    mFData <- funData::multiFunData(lapply(FPCA, function (mark) {
      funData::funData(argvals = mark$workGrid, X = fitted(mark))
    }))
    uniExpansions <- lapply(FPCA, function (mark) {
      list(type = "given", functions = funData::funData(
        argvals = mark$workGrid, X = t(mark$phi)))
    })
    if (is.null(M)) {
      M <- sum(sapply(FPCA, "[[", "selectK"))
    }

    # Extract weights
    weight_vec <- 1/sapply(lapply(FPCA, "[[", "lambda"), sum)

  } else if (method == "PACE") {

    # Construct irregular FunData
    m_irregFunData <- lapply(marker_dat, function (mark) {
      mark <- mark[order(mark[, time]), ]
      funData::irregFunData(argvals = split(mark[, time], mark[, id]),
                            X = split(mark$res, mark[, id]))
    })

    # Remove observations with too few scalar observations
    if (!is.null(npc)) {
      rem <- lapply(m_irregFunData, function (mark) {
        which(lapply(mark@argvals, length) < npc)
      })
      rem <- Reduce(union, rem)
      take <- if (length(rem) == 0) {
        seq_len(funData::nObs(m_irregFunData[[1]]))
      } else {
        seq_len(funData::nObs(m_irregFunData[[1]]))[-rem]
      }
      m_irregFunData <- lapply(m_irregFunData, function (mark) {
        funData::extractObs(mark, obs = take)
      })
    }

    # Use PACE function for each marker
    FPCA <- lapply(m_irregFunData, function(mark) {
      MFPCA::PACE(mark, npc = npc, pve = pve_uni, nbasis = nbasis)
    })

    # Construct multivariate FunData and estimated FPCs
    mFData <- funData::multiFunData(lapply(FPCA, "[[", "fit"))
    uniExpansions <- lapply(FPCA, function (mark) {
      list(type = "given", functions = mark$functions)
    })
    if (is.null(M)) {
      M <- sum(sapply(FPCA, "[[", "npc"))
    }

    # Extract weights
    weight_vec <- 1/sapply(lapply(FPCA, "[[", "values"), sum)
  }

  if (!weights) {
    weight_vec <- rep(1, length(mFData))
  }

  MFPCA <- MFPCA::MFPCA(mFData = mFData, M = M, uniExpansions = uniExpansions,
                        weights = weight_vec, fit = fit)
  attr(MFPCA, "sigma2") <- lapply(FPCA, "[[", "sigma2")
  if (save_uniFPCA) {
    attr(MFPCA, "uniFPCA") <- FPCA
  }
  if (save_uniGAM) {
    attr(MFPCA, "uniGAM") <- marker_mod
  }
  if (fit) {
    attr(MFPCA, "mFData") <- mFData
  }
  MFPCA
}


# Attach Weighted Functional Principal Components to the Data -----------


#' Attach Weighted Functional Principal Components to the Data
#'
#' @param mfpca MFPCA object from which to extract the weighted FPCS.
#' @param data Data set to which the weighted FPCS are to be attached.
#' @param n Number of FPCs to attach. Defaults to NULL which corresponds to all
#'  FPCs in mfpc.
#' @param obstime Name of the time variable in data set at which points to
#'  evaluate.
#' @param marker Name of the marker variable in the data set which separates the
#'  data.
#' @param eval_weight Weight the FPC by the square root of its eigenvalue (then
#'  variance comparable throughout all FPCs). Defaults to FALSE.
#' @export
#' @returns Data set supplied as argument \code{data} with additional columns
#'  corresponding to the evaluations of the MFPC basis.
#' @export
#' @examples
#' # Small example based on subset of PBC data
#' data(pbc_subset)
#'
#' # Estimate MFPC basis and attach to data
#' mfpca <- preproc_MFPCA(pbc_subset, uni_mean = paste0(
#'     "logy ~ 1 + sex + drug + s(obstime, k = 5, bs = 'ps') + ",
#'     "s(age, k = 5, bs = 'ps')"),
#'     pve_uni = 0.99, nbasis = 5, weights = TRUE, save_uniFPCA = TRUE)
#' pbc_subset <- attach_wfpc(mfpca, pbc_subset, n = 2)
attach_wfpc <- function(mfpca, data, n = NULL, obstime = "obstime",
                        marker = "marker", eval_weight = FALSE){

  # Is the data sorted by marker
  if (!all(order(data[[marker]]) == seq_len(nrow(data)))) message("ORDER!")

  wfpc <- NULL

  splitdat <- split(data[[obstime]], data[[marker]])

  if (!is.null(n)) {
    mfpca <- list(values = mfpca$values[seq_len(n)],
                  functions = funData::extractObs(mfpca$functions,
                                                  obs = seq_len(n)))
  }

  # eval_mpfc evaluates on all markers, so choose only the current one
  for (mark in seq_along(levels(data[[marker]]))) {
    mobs <- length(splitdat[[mark]])
    tot_wfpc <- eval_mfpc(mfpca = mfpca, timepoints = splitdat[[mark]],
                          eval_weight = eval_weight)
    wfpc <- rbind(wfpc, tot_wfpc[(mark-1)*mobs + seq_len(mobs), , drop = FALSE])
  }

  data <- cbind(data, wfpc)
  data

}
