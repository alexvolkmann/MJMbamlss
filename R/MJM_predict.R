## Prediction.
#' Prediction of MJM model
#'
#' Note: Writing a predict function is a bit tricky.
#' For longitudinal prediction, if subject specific predictions are wanted, then
#' the PCRE terms must be attached to newdata and already evaluated.
#' If the model uses standardized survival matrices, the different linear
#' predictors should be predicted using different data sets.
#' @param object bamlss-model object to be predicted.
#' @param newdata Dataset for which to create predictions. Not needed for
#'   conditional survival probabilities.
#' @param type Character string indicating which type of predictions to compute.
#'   \code{link} returns estimates for all predictors with the respective link
#'   functions applied, \code{"parameter"} returns the estimates for all
#'   pedictors, \code{"probabilities"} returns the survival probabilities
#'   conditional on the survival up to the last longitudinal measurement, and
#'   \code{"cumhaz"} return the cumulative hazard up to the survival time or for
#'   a time window after the last longitudinal measurement. If \code{type} is
#'   set to \code{"loglik"}, the log-likelihood of the joint model is returned.
#'   Note that types \code{"probabilities"} and \code{"cumhaz"} are not yet
#'   implemented.
#' @param id Integer or character, that specifies the individual for which the
#'   plot should be created.
#' @param dt The time window after the last observed measurement for which
#'   predictions should be computed.
#' @param FUN A function that should be applied on the samples of predictors or
#'   parameters, depending on argument \code{type}.
#' @param subdivisions Number of Gaussian quadrature points for survival
#'   integral calculation.
#' @param cores Specifies the number of cores that should be used for
#'   prediction. Note that this functionality is based on the
#'   \code{\link[parallel]{parallel}} package.
#' @param chunks Should computations be split into \code{chunks}? Prediction is
#'   then processed sequentially.
#' @param verbose Print information during runtime of the algorithm.
#' @param ... Currently not used.
MJM_predict <- function(object, newdata,
                        type = c("link", "parameter", "probabilities",
                                 "cumhaz"),
                       dt, id,
                       FUN = function(x) { mean(x, na.rm = TRUE) },
                       subdivisions = 7, cores = NULL,
                       chunks = 1, verbose = FALSE,  ...)
{

  if(missing(dt)) dt <- 0
  if(missing(id)) i <- NULL else i <- id
  idvar <- attr(object$y, "idvar")
  marker_name <- attr(object$y, "marker_name")
  nmarker <- attr(object$y, "nmarker")
  marker_levels <- levels(attr(object$y, "marker"))
  timevar_mu <- attr(object$y, "timevar")["mu"]


  if(length(type) > 1)
    type <- type[1]
  type <- match.arg(type)

  if(type == "probabilities"){
    if(!is.null(newdata)){
      warning("The provided newdata will be ignored for the prediction",
              "of conditional survival probabilities.")
      newdata <- NULL
    }
    if(dt == 0){
      stop("Please specify a time window for the prediction of conditional",
           "survival probabilities.")
    }
  }

  if(type == "cumhaz" & !is.null(newdata) & dt > 0){
    warning("The provided newdata will be ignored for the prediction of",
            "conditional survival t + dt.")
    newdata <- NULL
  }


  if(is.null(newdata)){
    newdata <- model.frame(object)
  }


  if(!(type %in% c("probabilities", "cumhaz"))) {

    object$family$predict <- NULL

    # Check for PCRE-term
    pcre_mu <- lapply(object$x$mu$smooth.construct, function (x) {
      any(grepl("pcre.random.effect", class(x), fixed = TRUE))
    })
    for (j in which(unlist(pcre_mu))) {

      # Compute the evaluated fpc basis functions for the newdata
      attr(object$x$mu$smooth.construct[[j]], "which_marker") <-
        newdata[, attr(object$y, "marker_name")]

      # Check and attach missing FPCbasis (set to 0)
      add_fpc <- which(!object$x$mu$smooth.construct[[j]]$term %in%
                         colnames(newdata))
      if (length(add_fpc) > 0) {
        for (add in add_fpc) {
          newdata[[object$x$mu$smooth.construct[[j]]$term[add]]] <- 0
        }
      }

    }

    return(bamlss:::predict.bamlss(object, newdata = newdata, type = type,
                                   FUN = FUN, cores = cores, chunks = chunks,
                                   verbose = verbose, ...))
  }

  stop("Types 'probabilities' and 'cumhaz' are not yet implemented.")
  # if(object$family$family != "mjm")
  #   stop("object must be a mjm-type model!")
  #
  # timevar <- attr(object$y, "timevar")
  # tmax_model <- max(newdata[,timevar["lambda"]])
  #
  # if(!is.null(i)){
  #   if(!is.character(i))
  #     i <- levels(newdata[[idvar]])[i]
  #   newdata <- subset(newdata, newdata[[idvar]] %in% i)
  # }
  #
  # tmax_pred <- max(newdata[,timevar["mu"]]) + dt
  #
  # if(tmax_pred > tmax_model){
  #   warning("Predictions should not be made beyond the modelled time range.
  #           Please adjust the time window dt accordingly.")
  # }
  #
  #
  #
  # ## Create the time grid.
  # # Gaussian Quadrature
  # stopifnot(requireNamespace("statmod"))
  # gq <- statmod::gauss.quad(subdivisions)
  #
  # grid <- function(lower, upper) {
  #   (upper - lower) / 2 * gq$nodes + (upper + lower) / 2
  # }
  #
  # jm_probs <- function(data) {
  #
  #   if(dt == 0){
  #     take <- !duplicated(data[, c(timevar["lambda"], idvar)])
  #     dsurv <- subset(data, take)
  #     timegrid <- lapply(dsurv[[timevar["lambda"]]],
  #                        function(x){grid(0, x)})
  #   } else {
  #     take <- !duplicated(data[, c(timevar["lambda"], idvar)], fromLast = TRUE)
  #     dsurv <- subset(data, take)
  #     timegrid <- lapply(dsurv[[timevar["mu"]]],
  #                        function(x){grid(x, x+dt)})
  #   }
  #   nobs <- nrow(dsurv)
  #   gdim <- c(length(timegrid), length(timegrid[[1]]))
  #
  #   # long data grid for multiple markers
  #   dsurv_long <- dsurv[rep(seq_len(nrow(dsurv)), nmarker), ]
  #   dsurv_long[[marker_name]] <- rep(marker_levels, each = nrow(dsurv))
  #   timegrid_long <- rep(timegrid, nmarker)
  #
  #   pred.setup <- bamlss:::predict.bamlss(object, data, type = "link",
  #                                         get.bamlss.predict.setup = TRUE, ...)
  #
  #   enames <- pred.setup$enames
  #
  #   pred_gamma <- with(pred.setup,
  #                      bamlss:::.predict.bamlss(
  #                        "gamma", object$x$gamma, samps, enames$gamma,
  #                        intercept, nsamps, dsurv))
  #
  #   pred_lambda <- with(pred.setup,
  #                       .predict.bamlss.mjm.td(
  #                         id = "lambda",
  #                         x = object$x$lambda$smooth.construct, samps = samps,
  #                         enames = enames$lambda, intercept = intercept,
  #                         nsamps = nsamps, newdata = dsurv,
  #                         yname = timevar["lambda"], grid = timegrid,
  #                         newdata_long = NULL, grid_long = NULL,
  #                         formula = bamlss:::drop.terms.bamlss(
  #                           object$x$lambda$terms, sterms = FALSE,
  #                           keep.response = FALSE), idvar = idvar,
  #                         timevar_mu = timevar_mu, nmarker = nmarker))
  #
  #   pred_mu <- with(pred.setup,
  #                   .predict.bamlss.mjm.td(
  #                     id = "mu", x = object$x$mu$smooth.construct,
  #                     samps = samps, enames = enames$mu, intercept = intercept,
  #                     nsamps = nsamps, newdata = dsurv,
  #                     yname = timevar["lambda"],
  #                     grid = timegrid,
  #                     newdata_long = dsurv_long, grid_long = timegrid_long,
  #                     formula = bamlss:::drop.terms.bamlss(
  #                       object$x$mu$terms, sterms = FALSE,
  #                       keep.response = FALSE), idvar = NULL,
  #                     timevar_mu = timevar_mu, nmarker = nmarker))
  #
  #   pred_alpha <- with(pred.setup,
  #                      .predict.bamlss.mjm.td(
  #                        id = "alpha", object$x$alpha$smooth.construct,
  #                        samps = samps, enames = enames$alpha,
  #                        intercept = intercept, nsamps = nsamps,
  #                        newdata = dsurv,  yname = timevar["lambda"],
  #                        grid = timegrid_long,
  #                        newdata_long = dsurv_long, grid_long = timegrid_long,
  #                        formula = bamlss:::drop.terms.bamlss(
  #                          object$x$alpha$terms, sterms = FALSE,
  #                          keep.response = FALSE), idvar = NULL,
  #                        timevar_mu = timevar_mu, nmarker = nmarker))
  #
  #   eta_timegrid_long <- colSums(aperm(array(pred_alpha * pred_mu,
  #                                            dim = c(nrow(pred_lambda), nmarker,
  #                                                    ncol(pred_lambda))),
  #                                      c(2, 1, 3)))
  #
  #   eta_timegrid <- pred_lambda + eta_timegrid_long
  #
  #   if(dt == 0){
  #
  #     # Only for cumulative hazard so integral starts at 0
  #     probs <- matrix(nrow = nobs,  ncol = ncol(eta_timegrid))
  #     for(i in seq_len(ncol(eta_timegrid))) {
  #       eta <- matrix(eta_timegrid[, i], nrow = gdim[1], ncol = gdim[2],
  #                     byrow = TRUE)
  #       eeta <- exp(eta)
  #
  #       # GQ integral
  #       int <- dsurv[[timevar["lambda"]]] / 2 * rowSums(t(gq$weights*t(eeta)))
  #
  #       probs[, i] <- exp(pred_gamma[, i]) * int
  #     }
  #
  #     if(!is.null(FUN)) {
  #       if(is.matrix(probs)) {
  #         if(ncol(probs) > 1)
  #           probs <- apply(probs, 1, FUN)
  #         probs <- t(probs)
  #       }
  #     }
  #
  #   } else {
  #     # Only starts after survival time
  #     probs <- matrix(nrow = nobs,  ncol = ncol(eta_timegrid))
  #     for(i in seq_len(ncol(eta_timegrid))) {
  #       eta <- matrix(eta_timegrid[, i], nrow = gdim[1], ncol = gdim[2],
  #                     byrow = TRUE)
  #       eeta <- exp(eta)
  #
  #       # GQ integral
  #       int <- dt / 2 * rowSums(t(gq$weights*t(eeta)))
  #
  #       probs[, i] <- exp(pred_gamma[, i]) * int
  #
  #     }
  #     if (type == "probabilities") {
  #       probs <- exp(-1 * probs)
  #     }
  #
  #     if(!is.null(FUN)) {
  #       if(is.matrix(probs)) {
  #         if(ncol(probs) > 1)
  #           probs <- apply(probs, 1, FUN)
  #         probs <- t(probs)
  #       }
  #     }
  #   }
  #
  #
  #   return(probs)
  # }
  #
  # ia <- interactive()
  #
  # if(is.null(cores)) {
  #   if(chunks < 2) {
  #     probs <- jm_probs(newdata)
  #   } else {
  #     id <- sort(rep(1:chunks, length.out = nrow(newdata)))
  #     newdata <- split(newdata, id)
  #     chunks <- length(newdata)
  #     probs <- NULL
  #     for(i in 1:chunks) {
  #       if(verbose) {
  #         cat(if(ia) "\r" else "\n")
  #         cat("predicting chunk", i, "of", chunks, "...")
  #         if(.Platform$OS.type != "unix" & ia) utils::flush.console()
  #       }
  #       if(i < 2) {
  #         probs <- jm_probs(newdata[[i]])
  #       } else {
  #         if(is.null(dim(probs))) {
  #           probs <- c(probs, jm_probs(newdata[[i]]))
  #         } else {
  #           probs <- rbind(probs, jm_probs(newdata[[i]]))
  #         }
  #       }
  #     }
  #     if(verbose) cat("\n")
  #   }
  # } else {
  #   parallel_fun <- function(i) {
  #     if(chunks < 2) {
  #       pr <- jm_probs(newdata[[i]])
  #     } else {
  #       idc <- sort(rep(1:chunks, length.out = nrow(newdata[[i]])))
  #       nd <- split(newdata[[i]], idc)
  #       chunks <- length(nd)
  #       pr <- NULL
  #       for(j in 1:chunks) {
  #         if(j < 2) {
  #           pr <- jm_probs(nd[[j]])
  #         } else {
  #           if(is.null(dim(pr))) {
  #             pr <- c(pr, jm_probs(nd[[j]]))
  #           } else {
  #             pr <- rbind(pr, jm_probs(nd[[j]]))
  #           }
  #         }
  #       }
  #     }
  #     return(pr)
  #   }
  #
  #   id <- sort(rep(1:cores, length.out = nrow(newdata)))
  #   newdata <- split(newdata, id)
  #   cores <- length(newdata)
  #   probs <- parallel::mclapply(1:cores, parallel_fun, mc.cores = cores)
  #
  #   probs <- if(is.matrix(probs[[1]])) {
  #     do.call("rbind", probs)
  #   } else {
  #     do.call("c", probs)
  #   }
  # }
  #
  # return(drop(probs))
}



# .predict.bamlss.mjm.td <- function(id, x, samps, enames, intercept, nsamps,
#                                    newdata, yname, grid, formula, idvar,
#                                    timevar_mu, nmarker, newdata_long,
#                                    grid_long)
# {
#
#   # id is "lambda", "mu" etc, so could be renamed
#   snames <- colnames(samps)
#   enames <- gsub("p.Intercept", "p.(Intercept)", enames, fixed = TRUE)
#   has_intercept <- any(grepl(paste(id, "p", "(Intercept)", sep = "."),
#                              snames, fixed = TRUE))
#
#   if(intercept & has_intercept)
#     enames <- c("p.(Intercept)", enames)
#   if (!has_intercept) {
#     enames <- enames[-grep("p.(Intercept)", enames, fixed = TRUE)]
#   }
#   enames <- unique(enames)
#   ec <- sapply(enames, function(x) {
#     paste(strsplit(x, "")[[1]][1:2], collapse = "")
#   })
#   enames2 <- sapply(enames, function(x) {
#     paste(strsplit(x, "")[[1]][-c(1:2)], collapse = "")
#   })
#
#   # Long time grid not necessary for lambda
#   if(id == "lambda") {
#     grid_long <- grid
#     newdata_long <- newdata
#   }
#
#   eta <- 0
#   p_components <- grep("p.", ec)
#   if(length(p_components)) {
#     intcpt <- ifelse(has_intercept, "1", "-1")
#     f <- as.formula(paste("~", intcpt, "+",
#                           paste(enames2[p_components], collapse = "+")))
#     X <- param_time_transform_mjm(x = list(), formula = f, data = newdata_long,
#                                   grid = grid_long, yname = yname,
#                                   timevar = if (id == "mu") {timevar_mu}
#                                   else {yname}, take = NULL,
#                                   idvar = idvar, y, timevar2 =
#                                     if (id == "lambda") {timevar_mu}
#                                   else {NULL})$Xgrid
#
#     sn <- snames[bamlss:::grep2(paste(id, "p",
#                                       enames2[p_components], sep = "."),
#                                 snames, fixed = TRUE)]
#     if(!length(sn))
#       sn <- snames[bamlss:::grep2(paste(id, "p.model.matrix",
#                                         enames2[p_components], sep = "."),
#                                   snames, fixed = TRUE)]
#     eta <- eta + bamlss:::fitted_matrix(X, samps[, sn, drop = FALSE])
#   }
#   if(length(i <- grep("s.", ec))) {
#     y2 <- newdata[, yname, drop = FALSE]
#     for(j in enames2[i]) {
#       for(jj in grep(j, names(x), fixed = TRUE, value = TRUE)) {
#
#         if(!inherits(x[[jj]], "no.mgcv") & !inherits(x[[jj]], "special")) {
#           if(inherits(x[[j]],
#                       "pcre.random.effect")) {
#             x[[jj]]$term <- x[[jj]]$term[-length(x[[jj]]$term)]
#             X <- sm_time_transform_mjm_pcre(
#               x = x[[jj]], data = newdata,
#               grid = grid, yname = yname, timevar = timevar_mu,
#               take = NULL, nmarker = nmarker)$Xgrid
#           } else {
#             X <- sm_time_transform_mjm(
#               x = x[[jj]], data = newdata_long,
#               grid = grid_long, yname = yname,
#               timevar = if (id == "mu") timevar_mu else yname,
#               take = NULL, y = y2)$Xgrid
#           }
#
#           sn <- snames[bamlss:::grep2(paste(id, "s", jj, sep = "."), snames,
#                              fixed = TRUE)]
#           random <- if(!is.null(x[[jj]]$margin)) {
#             any(sapply(x[[jj]]$margin, function(z) {
#               inherits(z, "random.effect")
#               }))
#           } else inherits(x[[jj]], "random.effect")
#           ok <- if(random) {
#             if(ncol(X) == ncol(samps[, sn, drop = FALSE])) TRUE else FALSE
#           } else TRUE
#           if(ok) {
#             eta <- eta + bamlss:::fitted_matrix(X, samps[, sn, drop = FALSE])
#           } else {
#             warning("predictions do not include ", jj)
#           }
#         } else {
#           stop("no predictions for special terms available yet!")
#         }
#       }
#     }
#   }
#   eta
# }
