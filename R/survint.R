#' Survival Integral
#'
#' This function is a wrapper function for calculating the survival integral in
#' C needed in the calculation of the score vector and Hessian.
#'
#' The survival integral has a similar structure for the different model
#' predictors. It is always a sum over all individuals, followed by the
#' multiplication with a pre-integral factor (pre_fac). For the gamma predictor
#' a pre-integral vector is next. Then, the integral itself consists of a
#' weighted sum (weights) of gauss-quadrature integration points weighted by
#' the survival time of the individuals (survtime). Inside the integral, the
#' current additive predictor (omega) is multiplied with an in-integral vector
#' (int_vec), except for predictor gamma. All longitudinal predictors
#' addtitionally include an in-integration factor (int_fac).
#'
#' The difference between predictors "long" and "fpc_re" is that the latter
#' makes efficient use of the block structure of the design matrix for
#' unconstrained functional principal component random effects. The outputs
#' also differ as the Hessian for "fpc_re" is a diagonal matrix, so only the
#' diagonal elements are returned.
#'
#' @param pred String to define for which predictor the survival integral is
#'   calculated.
#' @param pre_fac Vector serving as factor before the survival integral.
#'   Corresponds to the gamma predictor.
#' @param pre_vec Matrix serving as row vectors before the survival integral.
#'   Only needed if pred = "gamma".
#' @param omega Vector serving as additive predictor placeholder within the
#'   survival integral. Present for all pred.
#' @param int_fac Vector serving as factor within the survival integral. Only
#'   needed for the longitudinal predictors.
#' @param int_vec Matrix serving as row vectors within the survival integral.
#'   NULL only if pred = "gamma".
#' @param weights Vector containing the Gaussian integration weights.
#' @param survtime Vector containing the survival times for weighting of the
#'   integral.
#' @useDynLib MJMbamlss
survint_C <- function(pred = c("lambda", "gamma", "long", "fpc_re"),
  pre_fac, pre_vec = NULL, omega, int_fac = NULL,
  int_vec = NULL, weights, survtime)
{

  if (is.null(pre_vec)) pre_vec <- 0
  if (is.null(int_fac)) int_fac <- 0
  if (is.null(int_vec)) int_vec <- 0

  if (pred == "fpc_re") {
    .Call("survint_re", pre_fac, omega, int_fac, int_vec, weights,
          survtime)
  } else {
    pred <- switch(pred,
                   "lambda" = 1L,
                   "gamma" = 2L,
                   "long" = 3L
    )
    .Call("survint", pred, pre_fac, pre_vec, omega, int_fac, int_vec,
          weights, survtime)
  }

}

