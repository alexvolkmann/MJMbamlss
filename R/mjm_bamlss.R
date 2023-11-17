#' Flexible Joint Models for Multivariate Longitudinal and Time-to-Event Data
#'
#' This package contains all functions and implementations of the corresponding
#' paper by Volkmann, Umlauf, Greven: "Flexible joint models for multivariate
#' longitudinal and time-to-event data using multivariate functional principal
#' components". Code to reproduce the simulation and analysis as well as
#' additional information on the model fitting process are contained in
#' the "inst" folder.
#'
#' @docType package
#' @name varbinq
NULL

#' Family for Flexible Multivariate Joint Model
#'
#' This function specifies the different predictors and link functions as well
#' as the corresponding transform/updating/sampling functions as well as the
#' predict function.
#'
#' Family object to fit a flexible additive joint model for multivariate
#' longitudinal and survival data under a Bayesian approach using multivariate
#' functional principal components as presented in Volkmann, Umlauf, Greven
#' (2023).
#' @param ... All arguments are actually hard coded as needed by the
#'   implementation.
#' @export
#' @import bamlss
#' @import stats
#' @import coda
#' @import utils
#' @references Volkmann, A., Umlauf, N., Greven, S. (2023). Flexible joint
#'   models for multivariate longitudinal and time-to-event data using
#'   multivariate functional principal components. <arXiv:2311.06409>
#' @returns An object of class \code{family.bamlss}.
#' @examples
#' data(pbc_subset)
#' mfpca <- preproc_MFPCA(pbc_subset, uni_mean = paste0(
#'     "logy ~ 1 + sex + drug + s(obstime, k = 10, bs = 'ps') + ",
#'     "s(age, k = 10, bs = 'ps')"),
#'     pve_uni = 0.99, nbasis = 5, weights = TRUE, save_uniFPCA = TRUE)
#' pbc_subset <- attach_wfpc(mfpca, p_long, n = 3)
mjm_bamlss <- function(...)
{
  links = c(
    lambda = "log",
    gamma = "log",
    mu = "identity",
    sigma = "log",
    alpha = "identity"
  )

  rval <- list(
    "family" = "mjm",
    "names" = c("lambda", "gamma", "mu", "sigma", "alpha"),
    "links" = links,
    "transform" = MJM_transform,
    "optimizer" = MJM_opt,
    "sampler" = MJM_mcmc,
    "predict" = MJM_predict
  )

  class(rval) <- "family.bamlss"
  rval
}
