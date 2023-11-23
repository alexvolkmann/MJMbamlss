#' Functional principal components analysis by smoothed covariance
#'
#' Decomposes functional observations using functional principal components
#' analysis. A mixed model framework is used to estimate scores and obtain
#' variance estimates. This function is a slightly adapted copy of the
#' \code{\link[refund]{fpca.sc}} function in package \code{refund}
#' (version 0.1-30).
#'
#' This function computes a FPC decomposition for a set of observed curves,
#' which may be sparsely observed and/or measured with error. A mixed model
#' framework is used to estimate curve-specific scores and variances.
#'
#' FPCA via kernel smoothing of the covariance function, with the diagonal
#' treated separately, was proposed in Staniswalis and Lee (1998) and much
#' extended by Yao et al. (2005), who introduced the 'PACE' method.
#' \code{fpca.sc} uses penalized splines to smooth the covariance function, as
#' developed by Di et al. (2009) and Goldsmith et al. (2013). This
#' implementation uses REML and Cederbaum et al. (2018) for smoothing the
#' covariance function.
#'
#' The functional data must be supplied as either \itemize{ \item an \eqn{n
#' \times d} matrix \code{Y}, each row of which is one functional observation,
#' with missing values allowed; or \item a data frame \code{ydata}, with
#' columns \code{'.id'} (which curve the point belongs to, say \eqn{i}),
#' \code{'.index'} (function argument such as time point \eqn{t}), and
#' \code{'.value'} (observed function value \eqn{Y_i(t)}).}
#'
#' @param Y,ydata the user must supply either \code{Y}, a matrix of functions
#' observed on a regular grid, or a data frame \code{ydata} representing
#' irregularly observed functions. See Details.
#' @param Y.pred if desired, a matrix of functions to be approximated using
#' the FPC decomposition.
#' @param argvals the argument values of the function evaluations in \code{Y},
#'  defaults to a equidistant grid from 0 to 1.
#' @param argvals_obs Should the timepoints of the original observations be
#'    used to evaluate the FPCs. Defaults to FALSE.
#' @param argvals_pred Vector of timepoints on which to evaluate the FPCs.
#'   Defaults to a sequence from 0 to 1.
#' @param random.int If \code{TRUE}, the mean is estimated by
#' \code{\link[gamm4]{gamm4}} with random intercepts. If \code{FALSE} (the
#' default), the mean is estimated by \code{\link[mgcv]{gam}} treating all the
#' data as independent.
#' @param nbasis number of B-spline basis functions used for estimation of the
#' mean function.
#' @param nbasis_cov number of basis functions used for the bivariate
#' smoothing of the covariance surface.
#' @param bs_cov type of spline for the bivariate smoothing of the covariance
#' surface. Default is symmetric fast covariance smoothing proposed by
#' Cederbaum.
#' @param pve proportion of variance explained: used to choose the number of
#' principal components.
#' @param npc prespecified value for the number of principal components (if
#' given, this overrides \code{pve}).
#' @param useSymm logical, indicating whether to smooth only the upper
#' triangular part of the naive covariance (when \code{cov.est.method==2}).
#' This can save computation time for large data sets, and allows for
#' covariance surfaces that are very peaked on the diagonal.
#' @param makePD logical: should positive definiteness be enforced for the
#' covariance surface estimate?
#' @param center logical: should an estimated mean function be subtracted from
#' \code{Y}? Set to \code{FALSE} if you have already demeaned the data using
#' your favorite mean function estimate.
#' @param cov.est.method covariance estimation method. If set to \code{1}
#' (the default), a one-step method that applies a bivariate smooth to the
#' \eqn{y(s_1)y(s_2)} values. This can be very slow. If set to \code{2}, a two-step
#' method that obtains a naive covariance estimate which is then smoothed.
#' @param integration quadrature method for numerical integration; only
#' \code{'trapezoidal'} is currently supported.
#' @return An object of class \code{fpca} containing:
#' \item{Yhat}{FPC approximation (projection onto leading components)
#' of \code{Y.pred} if specified, or else of \code{Y}.}
#' \item{Y}{the observed data}\item{scores}{\eqn{n
#' \times npc} matrix of estimated FPC scores.} \item{mu}{estimated mean
#' function (or a vector of zeroes if \code{center==FALSE}).} \item{efunctions
#' }{\eqn{d \times npc} matrix of estimated eigenfunctions of the functional
#' covariance, i.e., the FPC basis functions.} \item{evalues}{estimated
#' eigenvalues of the covariance operator, i.e., variances of FPC scores.}
#' \item{npc }{number of FPCs: either the supplied \code{npc}, or the minimum
#' number of basis functions needed to explain proportion \code{pve} of the
#' variance in the observed curves.} \item{argvals}{argument values of
#' eigenfunction evaluations} \item{sigma2}{estimated measurement error
#' variance.} \item{diag.var}{diagonal elements of the covariance matrices for
#' each estimated curve.} \item{VarMats}{a list containing the estimated
#' covariance matrices for each curve in \code{Y}.} \item{crit.val}{estimated
#' critical values for constructing simultaneous confidence intervals.}
#' @author Jeff Goldsmith \email{jeff.goldsmith@@columbia.edu}, Sonja Greven
#' \email{sonja.greven@@stat.uni-muenchen.de}, Lan Huo
#' \email{Lan.Huo@@nyumc.org}, Lei Huang \email{huangracer@@gmail.com}, and
#' Philip Reiss \email{phil.reiss@@nyumc.org}, Alexander Volkmann
#' @references Cederbaum, J. Scheipl, F. and Greven, S. (2018). Fast symmetric
#' additive covariance smoothing. \emph{Computational Statistics & Data
#' Analysis}, 120, 25--41.
#'
#' Di, C., Crainiceanu, C., Caffo, B., and Punjabi, N. (2009).
#' Multilevel functional principal component analysis. \emph{Annals of Applied
#' Statistics}, 3, 458--488.
#'
#' Goldsmith, J., Greven, S., and Crainiceanu, C. (2013). Corrected confidence
#' bands for functional data using principal components. \emph{Biometrics},
#' 69(1), 41--51.
#'
#' Staniswalis, J. G., and Lee, J. J. (1998). Nonparametric regression
#' analysis of longitudinal data. \emph{Journal of the American Statistical
#' Association}, 93, 1403--1418.
#'
#' Yao, F., Mueller, H.-G., and Wang, J.-L. (2005). Functional data analysis
#' for sparse longitudinal data. \emph{Journal of the American Statistical
#' Association}, 100, 577--590.
fpca <- function(Y = NULL, ydata = NULL, Y.pred = NULL, argvals = NULL,
                 argvals_obs = FALSE, argvals_pred = seq(0, 1, by = 0.01),
                 random.int = FALSE, nbasis = 10, nbasis_cov = nbasis,
                 bs_cov = "symm", pve = 0.99, npc = NULL,
                 useSymm = FALSE, makePD = FALSE, center = TRUE,
                 cov.est.method = 1, integration = "trapezoidal") {

  stopifnot((!is.null(Y) && is.null(ydata)) || (is.null(Y) && !is.null(ydata)))

  # if data.frame version of ydata is provided
  sparseOrNongrid <- !is.null(ydata)
  if (sparseOrNongrid) {
    stopifnot(ncol(ydata) == 3)
    stopifnot(c(".id", ".index", ".value") == colnames(ydata))
    stopifnot(is.null(argvals))
    irreg2mat <- utils::getFromNamespace("irreg2mat", "refund")
    Y = irreg2mat(ydata)
    argvals = sort(unique(ydata$.index))
    if (!argvals_obs) {
      argvals_pred <- sort(unique(argvals_pred))
    } else {
      argvals_pred <- argvals
    }
  }

  if (is.null(Y.pred))
    Y.pred = Y
  D = NCOL(Y)
  I = NROW(Y)
  I.pred = NROW(Y.pred)

  if (is.null(argvals)) {
    argvals = seq(0, 1, length = D)
    if (!argvals_obs) {
      argvals_pred <- argvals
    }
  }

  D_pred <- length(argvals_pred)

  d.vec = rep(argvals, each = I)
  id = rep(1:I, rep(D, I))

  if (center) {
    if (random.int) {
      ri_data <- data.frame(y = as.vector(Y), d.vec = d.vec, id = factor(id))
      gam0 = gamm4(y ~ s(d.vec, k = nbasis), random = ~(1 | id),
                   data = ri_data)$gam
      rm(ri_data)
    } else gam0 = gam(as.vector(Y) ~ s(d.vec, k = nbasis))
    mu = predict(gam0, newdata = data.frame(d.vec = argvals))
    muhat <- predict(gam0, newdata = data.frame(d.vec = argvals_pred))
    Y.tilde = Y - matrix(mu, I, D, byrow = TRUE)
  } else {
    Y.tilde = Y
    mu = rep(0, D)
  }

  if (cov.est.method == 2) {
    # smooth raw covariance estimate
    cov.sum = cov.count = cov.mean = matrix(0, D, D)
    for (i in 1:I) {
      obs.points = which(!is.na(Y[i, ]))
      cov.count[obs.points, obs.points] = cov.count[obs.points, obs.points] +
        1
      cov.sum[obs.points, obs.points] = cov.sum[obs.points, obs.points] +
        tcrossprod(Y.tilde[i, obs.points])
    }
    G.0 = ifelse(cov.count == 0, NA, cov.sum/cov.count)
    diag.G0 = diag(G.0)
    diag(G.0) = NA
    if (!useSymm) {
      row.vec = rep(argvals, each = D)
      col.vec = rep(argvals, D)
      npc.0 = matrix(predict(gam(as.vector(G.0) ~ te(row.vec, col.vec,
                                                     k = nbasis),
                                 weights = as.vector(cov.count)),
                             newdata = data.frame(row.vec = row.vec,
                                                  col.vec = col.vec)), D, D)
      npc.0 = (npc.0 + t(npc.0))/2
    } else {
      use <- upper.tri(G.0, diag = TRUE)
      use[2, 1] <- use[ncol(G.0), ncol(G.0) - 1] <- TRUE
      usecov.count <- cov.count
      usecov.count[2, 1] <- usecov.count[ncol(G.0), ncol(G.0) - 1] <- 0
      usecov.count <- as.vector(usecov.count)[use]
      use <- as.vector(use)
      vG.0 <- as.vector(G.0)[use]
      row.vec <- rep(argvals, each = D)[use]
      col.vec <- rep(argvals, times = D)[use]
      mCov <- gam(vG.0 ~ te(row.vec, col.vec, k = nbasis),
                  weights = usecov.count)
      npc.0 <- matrix(NA, D, D)
      spred <- rep(argvals, each = D)[upper.tri(npc.0, diag = TRUE)]
      tpred <- rep(argvals, times = D)[upper.tri(npc.0, diag = TRUE)]
      smVCov <- predict(mCov, newdata = data.frame(row.vec = spred,
                                                   col.vec = tpred))
      npc.0[upper.tri(npc.0, diag = TRUE)] <- smVCov
      npc.0[lower.tri(npc.0)] <- t(npc.0)[lower.tri(npc.0)]
    }
  } else if (cov.est.method == 1) {
    # smooth y(s1)y(s2) values to obtain covariance estimate
    row.vec = col.vec = G.0.vec = c()
    cov.sum = cov.count = cov.mean = matrix(0, D_pred, D_pred)
    for (i in 1:I) {
      obs.points = which(!is.na(Y[i, ]))
      temp = tcrossprod(Y.tilde[i, obs.points])
      diag(temp) = NA
      row.vec = c(row.vec, rep(argvals[obs.points], each = length(obs.points)))
      col.vec = c(col.vec, rep(argvals[obs.points], length(obs.points)))
      G.0.vec = c(G.0.vec, as.vector(temp))
      # still need G.O raw to calculate to get the raw to get the diagonal
      cov.count[obs.points, obs.points] = cov.count[obs.points, obs.points] +
        1
      cov.sum[obs.points, obs.points] = cov.sum[obs.points, obs.points] +
        tcrossprod(Y.tilde[i, obs.points])
    }
    row.vec.pred = rep(argvals_pred, each = D_pred)
    col.vec.pred = rep(argvals_pred, D_pred)
    npc.0 = matrix(predict(mgcv::bam(G.0.vec ~ s(row.vec, col.vec,
                                                 k = nbasis_cov, bs = bs_cov),
                                     method = "fREML"),
                           newdata = data.frame(row.vec = row.vec.pred,
                                                col.vec = col.vec.pred)),
                   D_pred, D_pred)
    #npc.0 = (npc.0 + t(npc.0))/2
    G.0 = ifelse(cov.count == 0, NA, cov.sum/cov.count)
    diag.G0 = diag(G.0)
  }

  if (makePD) {
    npc.0 <- {
      tmp <- Matrix::nearPD(npc.0, corr = FALSE, keepDiag = FALSE,
                            do2eigen = TRUE, trace = TRUE)
      as.matrix(tmp$mat)
    }
  }
  ### numerical integration for calculation of eigenvalues (see Ramsay &
  ### Silverman, Chapter 8)
  quadWeights <- utils::getFromNamespace("quadWeights", "refund")
  w <- quadWeights(argvals_pred, method = integration)
  # Wsqrt <- diag(sqrt(w))
  # Winvsqrt <- diag(1/(sqrt(w)))
  V <- t(sqrt(w)*t(sqrt(w)*npc.0))
    # Wsqrt %*% npc.0 %*% Wsqrt
  eigenV <- eigen(V, symmetric = TRUE)
  evalues = eigenV$values
  ###
  evalues = replace(evalues, which(evalues <= 0), 0)
  npc = ifelse(is.null(npc), min(which(cumsum(evalues)/sum(evalues) > pve)),
               npc)
  efunctions = matrix(1/sqrt(w)*
                        eigenV$vectors[, seq(len = npc)],
                      nrow = D_pred, ncol = npc)
  evalues = eigenV$values[1:npc]
  # use correct matrix for eigenvalue problem
  cov.hat = efunctions %*% tcrossprod(diag(evalues, nrow = npc, ncol = npc),
                                      efunctions)
  ### numerical integration for estimation of sigma2
  T.len <- argvals_pred[D_pred] - argvals[1]
  # total interval length
  T1.min <- min(which(argvals >= argvals[1] + 0.25 * T.len))
  # left bound of narrower interval T1
  T1.max <- max(which(argvals <= argvals_pred[D_pred] - 0.25 * T.len))
  # right bound of narrower interval T1
  DIAG = (diag.G0 - diag(cov.hat))[T1.min:T1.max]
  # function values
  w2 <- quadWeights(argvals[T1.min:T1.max], method = integration)
  sigma2 <- max(weighted.mean(DIAG, w = w2, na.rm = TRUE), 0)

  ####
  D.inv = diag(1/evalues, nrow = npc, ncol = npc)
  Z = efunctions
  Y.tilde = Y.pred - matrix(mu, I.pred, D, byrow = TRUE)
  Yhat = matrix(0, nrow = I.pred, ncol = D_pred)
  rownames(Yhat) = rownames(Y.pred)
  colnames(Yhat) = argvals_pred
  scores = matrix(NA, nrow = I.pred, ncol = npc)
  #VarMats = vector("list", I.pred)
  #for (i in 1:I.pred) VarMats[[i]] = matrix(NA, nrow = D_pred, ncol = D_pred)
  diag.var = matrix(NA, nrow = I.pred, ncol = D_pred)
  crit.val = rep(0, I.pred)
  for (i.subj in 1:I.pred) {
    obs.points = which(!is.na(Y.pred[i.subj, ]))
    if (sigma2 == 0 & length(obs.points) < npc)
      stop("Measurement error estimated to be zero and there are fewer",
          " observed points than PCs; scores cannot be estimated.")
    Zcur = matrix(Z[obs.points, ], nrow = length(obs.points), ncol = dim(Z)[2])
    ZtZ_sD.inv = solve(crossprod(Zcur) + sigma2 * D.inv)
    scores[i.subj, ] = ZtZ_sD.inv %*% t(Zcur) %*% (Y.tilde[i.subj, obs.points])
    Yhat[i.subj, ] = t(as.matrix(muhat)) + scores[i.subj, ] %*% t(efunctions)
    # if (var) {
    #   VarMats[[i.subj]] = sigma2 * Z %*% ZtZ_sD.inv %*% t(Z)
    #   diag.var[i.subj, ] = diag(VarMats[[i.subj]])
    #   if (simul & sigma2 != 0) {
    #     norm.samp = mvrnorm(2500, mu = rep(0, D), Sigma = VarMats[[i.subj]])/
    #       matrix(sqrt(diag(VarMats[[i.subj]])), nrow = 2500, ncol = D,
    #              byrow = TRUE)
    #     crit.val[i.subj] = quantile(apply(abs(norm.samp), 1, max), sim.alpha)
    #   }
    # }
  }

  ret.objects = c("Yhat", "Y", "scores",
                  if(argvals_obs) "mu" else "muhat",
                  "efunctions", "evalues", "npc",
                  if(argvals_obs) "argvals" else "argvals_pred",
                  "sigma2")

  ret = lapply(1:length(ret.objects), function(u) get(ret.objects[u]))
  names(ret) = ret.objects
  class(ret) = "fpca"
  return(ret)
}
