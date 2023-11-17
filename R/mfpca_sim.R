#' Function to calculate the multivariate FPCA for a given covariance matrix and
#' univariate basis functions
#'
#' @param cov Covariance matrix of the basis functions coefficients.
#' @param basis_funs List with basis functions on each dimension. The basis
#'  functions are funData objects
#' @param scores Matrix (n rows, B columns) containing the basis functions
#'  coefficients. Defaults to NULL which does not calculate the multivariate
#'  scores.
#' @param weights Vector of weights, defaults to 1 for each element
#' @export
#' @importFrom Matrix bdiag
#' @import foreach
#' @returns List mimicking an \code{MFPCAfit} object containing the following
#'   components:
#'   \describe{
#'    \item{values}{A vector of eigenvalues.}
#'    \item{functions}{A \code{multiFunData} object containing the multivariate
#'      functional principal components.}
#'    \item{scores}{A matrix containing the scores (if applicable).}
#'    \item{vectors}{A matrix representing the eigenvectors associated with the
#'      combined univaraite score vectors.}
#'    \item{normFactors}{The normalizing factors used for calculating the
#'      multivariate eigenfunctions and scores.}
#' }
#' @examples
#' library(funData)
#' # Covariance matrix for the data generation in simulation scenario I
#' auto <- matrix(c(0.08, -0.07, -0.07, 0.9), ncol = 2)
#' cross <- matrix(rep(0.03, 4), ncol = 2)
#' cor <- matrix(c(0, 1, 0.75, 0.5, 0, 0,
#'                 1, 0, 1, 0.75, 0.5, 0,
#'                 0.75, 1, 0, 1, 0.75, 0.5,
#'                 0.5, 0.75, 1, 0, 1, 0.75,
#'                 0, 0.5, 0.75, 1, 0, 1,
#'                 0, 0, 0.5, 0.75, 1, 0),
#'                ncol = 6)
#' cov <- kronecker(cor, cross) + kronecker(diag(c(1, 1.2, 1.4, 1.6, 1.8, 2)),
#'                                          auto)
#' # Basis functions on each dimension
#' seq1 <- seq(0, 1, by = 0.01)
#' b_funs <- rep(list(funData(argvals = seq1,
#'               X = matrix(c(rep(1, length(seq1)), seq1),
#'               byrow = TRUE, ncol = length(seq1)))), 6)
#' # Prepare objects for the model on different data sets
#' mfpca_tru <- MFPCA_cov(cov = cov, basis_funs = b_funs)
MFPCA_cov <- function (cov, basis_funs, scores = NULL, weights = NULL) {

  # To not get a CMD check note for undefined global function
  j <- NULL

  # Information about the objects
  p <- length(basis_funs)
  npc <- sapply(basis_funs, funData::nObs)
  npcCum <- cumsum(c(0, npc))
  argvals_list <- lapply(basis_funs, function(x) {
    argvals(x)[[1]]
  })
  if (is.null(weights)) {
    weights <- rep(1, p)
  }
  allWeights <- rep(sqrt(weights), npc)

  # Construct matrix of basis function integrals
  calcBasisIntegrals <- utils::getFromNamespace("calcBasisIntegrals", "MFPCA")
  B_block <- Matrix::bdiag(lapply(basis_funs, function (x){
    calcBasisIntegrals(x@X, dimSupp = 1, x@argvals)
  }))
  # Construct weighted coefficient block matrix
  Qw <- cov * outer(allWeights, allWeights, "*")

  # Eigendecomposition
  eig <- eigen(B_block %*% Qw)
  values <- Re(eig$values)
  vectors <- Re(eig$vectors)

  # before the sums
  normFactors <- 1/sqrt(diag(t(vectors) %*% Qw %*% vectors))
  # after the sums
  blockWeights <- Qw %*% vectors

  # Calculation of MFPCs
  eFunctions <- funData::multiFunData(
    foreach::'%do%'(foreach::foreach(j = seq_len(p)), {
      MFPCA::univExpansion(type = "given",
                            scores = 1/sqrt(weights[j] * values) * normFactors *
                              t(blockWeights[npcCum[j] + seq_len(npc[j]), ,
                                             drop = FALSE]),
                            argvals = argvals_list[j],
                            functions = basis_funs[[j]],
                            params = NULL)
  }))

  # Calculate the multivariate scores
  if (!is.null(scores)) {
    m_scores <- as.matrix((scores * sqrt(allWeights)) %*% vectors %*%
                            diag(sqrt(values) * normFactors))
  }

  out <- list("values" = values,
              "functions" = eFunctions,
              "scores" = if (is.null(scores)) NULL else m_scores,
              "vectors" = vectors,
              "normFactors" = normFactors)
}
