#' PBC Subset
#'
#' A subset of the pbc data provided by package \code{JMbayes2} used only to
#' illustrate the functions.
#'
#' @format ## `pbc_subset`
#' A data frame with 336 observations and 10 columns:
#' \describe{
#'    \item{id}{Subject id.}
#'    \item{survtime}{Survival time of composite endpoint.}
#'    \item{event}{Composite endpoint.}
#'    \item{sex}{Male or female.}
#'    \item{drug}{Placebo or D-penicil.}
#'    \item{age}{Age.}
#'    \item{marker}{Name of longitudinal biomarker (albumin, SerBilir, serChol,
#'      SGOT)}
#'    \item{obstime}{Longitudinal time.}
#'    \item{y}{Longitudinal outcome value.}
#'    \item{logy}{Log-transformed longitudinal outcome value.}
#' }
#' @source <https://cran.r-project.org/web/packages/JMbayes2/index.html>
"pbc_subset"
