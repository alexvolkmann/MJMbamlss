
# MJM_transform -----------------------------------------------------------

MJM_transform <- function(object, subdivisions = 7, timevar = NULL, tau = NULL,
                          idvar = NULL, uni = FALSE, std_surv = TRUE, ...) {


  # Gaussian Quadrature
  stopifnot(requireNamespace("statmod"))
  gq <- statmod::gauss.quad(subdivisions)


  # Get idvar
  ## The last term in mu has to be a random effect with the idvar in the first
  ## place
  class_mu <- lapply(object$x$mu$smooth.construct, class)
  class_mu <- sapply(class_mu, function(x) x[1])
  if(!any(class_mu %in% c("pcre.random.effect", "unc_pcre.random.effect"))){
    j <- length(class_mu)
  } else {
    j <- which(class_mu %in% c("pcre.random.effect",
                               "unc_pcre.random.effect"))[1]
  }
  if (is.null(idvar)) {
    smj <- object$x$mu$smooth.construct[[j]]
    idvar <- smj$term[1]
    message("ID taken to be ", idvar)
  }

  # Get marker variable (fixed for now) - for MJM_predict()
  marker_name <- "marker"
  if(uni == TRUE) {
    object$model.frame[, marker_name] <- factor(1)
  }

  # Get ordering of the design matrix by id for PCRE crossproduct
  # and the number of observations per individual
  crossprod_info <- list(order(object$model.frame[, idvar]),
                         as.integer(table(object$model.frame[, idvar])))

  # Function to transform intervals to [0, 1] for GQ integration
  create_grid <- function(time) {
    time / 2 * gq$nodes + time / 2
  }


  # Different helpful variables and grid of timepoints for GQ
  y2_l <- cbind(object$y[[1]][, "time"], object$model.frame[[idvar]])
  colnames(y2_l) <- c("time", idvar)
  take <- !duplicated(y2_l[, 1:2])
  take_last <- !duplicated(y2_l, fromLast = TRUE)
  nsubj <- length(unique(y2_l[, idvar]))
  y2 <- y2_l[take, , drop = FALSE]
  grid <- lapply(y2[, "time"], create_grid)

  yname <- all.names(object$x$lambda$formula[2])[2]
  timevar_mu <- timevar
  if(is.null(timevar_mu))
    stop("the time variable is not specified, needed for mu!")
  timevar <- yname

  if (!identical(seq_len(nrow(object$model.frame)),
      order(object$model.frame[, marker_name],
            object$model.frame[, idvar],
            object$model.frame[, timevar_mu]))) {
    stop("Order the data by marker-id-obstime!")
  }

  # Longitudinal variables
  nmarker <- 1
  marker <- rep(1, nrow(y2_l))
  if(!is.null(object$model.frame$marker)) {
    y2_l <- cbind(y2_l, "marker" = as.factor(object$model.frame$marker))
    nmarker <- length(levels(as.factor(object$model.frame$marker)))
    marker <- object$model.frame$marker
  }
  take_l <- !duplicated(y2_l)
  take_last_l <- !duplicated(y2_l, fromLast = TRUE)
  y2_l <- y2_l[take_l, , drop = FALSE]
  grid_l <- lapply(y2_l[, "time"], create_grid)


  # Save information for optimizer in attributes of y
  attr(object$y, "gq_weights") <- gq$weights
  attr(object$y, "status") <- object$y[[1]][, "status"][take_last]
  attr(object$y, "take_last") <- take_last
  attr(object$y, "take_last_l") <- take_last_l
  attr(object$y, "nsubj") <- nsubj
  attr(object$y, "nmarker") <- nmarker
  attr(object$y, "marker") <- marker
  attr(object$y, "marker_name") <- marker_name
  attr(object$y, "timevar") <- c("lambda" = timevar, "mu" = timevar_mu)
  attr(object$y, "idvar") <- idvar


  # DESIGN CONSTRUCT
  # design.construct als Funktion eingebaut, damit sich die eta-Vektoren
  # der Survival-Prädiktoren in der Länge ändern
  ## Recompute design matrixes for lambda, gamma, alpha.
  for(j in c("lambda", "gamma", "alpha")) {
    which_take <- if (j == "alpha") take_l else take
    object$x[[j]] <- bamlss:::design.construct(
      object$terms,
      data = object$model.frame[which_take, , drop = FALSE],
      model.matrix = TRUE,
      smooth.construct = TRUE,
      model = j,
      scale.x = FALSE)[[j]]
  }

  ## Add small constant to penalty diagonal for stabilization
  for(j in names(object$x)) {
    for(sj in names(object$x[[j]]$smooth.construct)) {
      if(length(object$x[[j]]$smooth.construct[[sj]]$S)){
        for(k in seq_along(object$x[[j]]$smooth.construct[[sj]]$S)){
          if(!is.list(object$x[[j]]$smooth.construct[[sj]]$S[[k]]) &
             !is.function(object$x[[j]]$smooth.construct[[sj]]$S[[k]])) {
            nc <- ncol(object$x[[j]]$smooth.construct[[sj]]$S[[k]])
            object$x[[j]]$smooth.construct[[sj]]$S[[k]] <-
              object$x[[j]]$smooth.construct[[sj]]$S[[k]] + diag(1e-08, nc, nc)
          }
        }
      }
    }
  }

  ## The basic setup.
  if(is.null(attr(object$x, "bamlss.engine.setup")))
    object$x <- bamlss.engine.setup(object$x, ...)


  ## Remove intercept from lambda.
  if(!is.null(object$x$lambda$smooth.construct$model.matrix)) {
    cn <- colnames(object$x$lambda$smooth.construct$model.matrix$X)
    if("(Intercept)" %in% cn) {
      object$x$lambda$smooth.construct$model.matrix$X <-
        object$x$lambda$smooth.construct$model.matrix$X[, cn != "(Intercept)",
                                                        drop = FALSE]
    }
    if(ncol(object$x$lambda$smooth.construct$model.matrix$X) < 1) {
      object$x$lambda$smooth.construct$model.matrix <- NULL
      object$x$lambda$terms <- bamlss:::drop.terms.bamlss(
        object$x$lambda$terms, pterms = FALSE, keep.intercept = FALSE)
    } else {
      object$x$lambda$smooth.construct$model.matrix$term <-
        gsub("(Intercept)+", "",
             object$x$lambda$smooth.construct$model.matrix$term, fixed = TRUE)
      object$x$lambda$smooth.construct$model.matrix$state$parameters <-
        object$x$lambda$smooth.construct$model.matrix$state$parameters[-1]
      object$x$lambda$terms <- bamlss:::drop.terms.bamlss(
        object$x$lambda$terms, pterms = TRUE, sterms = TRUE,
        keep.intercept = FALSE)
    }
  }


  ## Compute new design matrices for integration for smooth terms.
  for(i in c("lambda", "mu", "alpha")) {
    if(!is.null(object$x[[i]]$smooth.construct)) {
      for(j in names(object$x[[i]]$smooth.construct)) {
        if(j != "model.matrix") {
          xterm <- object$x[[i]]$smooth.construct[[j]]$term
          by <- if(object$x[[i]]$smooth.construct[[j]]$by != "NA") {
            object$x[[i]]$smooth.construct[[j]]$by
          } else NULL
          if(inherits(object$x[[i]]$smooth.construct[[j]],
                      "pcre.random.effect")) {
            unique_term <- unique(c(xterm, yname, by, timevar_mu, idvar))
            object$x[[i]]$smooth.construct[[j]] <-
              sm_time_transform_mjm_pcre(
                object$x[[i]]$smooth.construct[[j]],
                object$model.frame[, unique_term, drop = FALSE],
                grid, yname, timevar_mu, take_last, nmarker = nmarker,
                cp_info = crossprod_info)
          } else {
            unique_term <- unique(c(xterm, yname, by,
                                    if(i == "mu") timevar_mu else timevar,
                                    idvar))
            object$x[[i]]$smooth.construct[[j]] <- sm_time_transform_mjm(
              x = object$x[[i]]$smooth.construct[[j]],
              data = object$model.frame[, unique_term, drop = FALSE],
              grid = if(i == "lambda") grid else grid_l,
              yname = yname,
              timevar = if(i == "mu") timevar_mu else timevar,
              take = if(i == "lambda") take_last else take_last_l,
              y = if (i == "lambda") y2 else y2_l)
          }
        }
      }
    }
  }

  # Standardize the survival covariates
  if (std_surv) {
    object$x$gamma$smooth.construct$model.matrix$w_bar <-
      colMeans(object$x$gamma$smooth.construct$model.matrix$X)[-1]
    object$x$gamma$smooth.construct$model.matrix$w_sd <-
      apply(object$x$gamma$smooth.construct$model.matrix$X, 2, sd)[-1]
    object$x$gamma$smooth.construct$model.matrix$X[, -1] <- scale(
      object$x$gamma$smooth.construct$model.matrix$X[, -1],
      center = object$x$gamma$smooth.construct$model.matrix$w_bar,
      scale = object$x$gamma$smooth.construct$model.matrix$w_sd
    )
    object$std_surv <- TRUE
  } else {
    object$std_surv <- NULL
  }

  ## Now linear part.
  for(i in c("lambda", "mu", "alpha")) {
    if(!is.null(object$x[[i]]$smooth.construct$model.matrix)) {
      object$x[[i]]$smooth.construct$model.matrix <- param_time_transform_mjm(
        object$x[[i]]$smooth.construct$model.matrix,
        bamlss:::drop.terms.bamlss(object$x[[i]]$terms, sterms = FALSE,
                                   keep.response = FALSE), object$model.frame,
        if(i == "lambda") grid else grid_l, yname,
        if(i != "mu") timevar else timevar_mu,
        if(i == "lambda") take_last else take_last_l,
        idvar = if (i == "lambda") idvar else NULL,
        timevar2 = if (i == "lambda") timevar_mu else NULL)
    }
  }

  ## Update prior/grad/hess functions
  for(j in names(object$x)) {
    for(sj in names(object$x[[j]]$smooth.construct)) {
      if (!is.null(tau[[j]][[sj]])) {
        tau_pos <- grepl(
          "tau", names(object$x[[j]]$smooth.construct[[sj]]$state$parameters))
        if (sum(tau_pos) != length(tau[[j]][[sj]])) {
          warning(cat("Smooth", sj, "in", j, "does not have correct length.\n"))
        }
        object$x[[j]]$smooth.construct[[sj]]$state$parameters[tau_pos] <-
          tau[[j]][[sj]]
      }
      priors <- bamlss:::make.prior(object$x[[j]]$smooth.construct[[sj]])
      object$x[[j]]$smooth.construct[[sj]]$prior <- priors$prior
      object$x[[j]]$smooth.construct[[sj]]$grad <- priors$grad
      object$x[[j]]$smooth.construct[[sj]]$hess <- priors$hess
    }
  }


  # Include a function to compute the logLikelihood for a set of parameters
  object$family$p2logLik <- function (par, logPost = FALSE, ...) {
    eta <- eta_timegrid <- list()
    lprior <- 0
    nx <- c("lambda", "mu", "gamma", "sigma", "alpha")
    for(j in nx) {
      eta[[j]] <- 0
      if(!(j %in% c("gamma", "sigma"))) {
        eta_timegrid[[j]] <- 0
        if (j == "mu") {
          eta_T_mu <- 0
        }
      }

      for(sj in names(object$x[[j]]$smooth.construct)) {
        pn <- paste(j, if(sj != "model.matrix") "s" else "p", sep = ".")
        cn <- colnames(object$x[[j]]$smooth.construct[[sj]]$X)
        if(is.null(cn))
          cn <- paste("b", 1:ncol(object$x[[j]]$smooth.construct[[sj]]$X),
                      sep = "")
        pn0 <- paste(pn, sj, sep = ".")
        pn <- paste(pn0, cn, sep = ".")
        if(all(is.na(par[pn]))) {
          if(sj == "model.matrix")
            pn <- gsub(".p.model.matrix.", ".p.", pn, fixed = TRUE)
        }
        eta[[j]] <- eta[[j]] +
          object$x[[j]]$smooth.construct[[sj]]$X %*% par[pn]
        if(!(j %in% c("gamma", "sigma"))) {
          eta_timegrid[[j]] <- eta_timegrid[[j]] +
            object$x[[j]]$smooth.construct[[sj]]$Xgrid %*% par[pn]
          if (j == "mu") {
            eta_T_mu + object$x[[j]]$smooth.construct[[sj]]$XT %*% par[pn]
          }
        }

        if(logPost) {
          pn2 <- paste(pn0, "tau2", sep = ".")
          tpar <- par[bamlss:::grep2(c(pn, pn2), names(par), fixed = TRUE)]
          lprior <- lprior + object$x[[j]]$smooth.construct[[sj]]$prior(tpar)
        }
      }
    }
    eta_timegrid_long <- drop(
      t(rep(1, nmarker)) %x% diag(length(eta_timegrid$lambda)) %*%
        (drop(eta_timegrid$alpha) * drop(eta_timegrid$mu)))
    eta_timegrid <- eta_timegrid$lambda + eta_timegrid_long

    eta_T_long <- drop(
      t(rep(1, nmarker)) %x% diag(nsubj) %*% (eta$alpha*eta_T_mu))
    eta_T <- eta$lambda + eta$gamma + eta_T_long

    sum_Lambda <- drop(y2[, 1]/2 * exp(eta$gamma)) %*%
      (diag(nsubj)%x%t(gq$weights))%*%
      exp(eta_timegrid)
    logLik <- drop(object$y[[1]][, "status"][take_last] %*% eta_T- sum_Lambda) +
      sum(dnorm(object$y[[1]][, "obs"], mean = eta$mu, sd = exp(eta$sigma),
                log = TRUE))
    if(logPost)
      logLik <- logLik + lprior
    return(drop(logLik))
  }

  return(object)
}



# Smooth time transformer function ----------------------------------------

#' @importFrom mgcv PredictMat
sm_time_transform_mjm <- function(x, data, grid, yname, timevar, take, y) {

  if(!is.null(take))
    data <- data[take, , drop = FALSE]
  X <- XT <- NULL
  for(j in x$term) {
    if((j != yname) & (j != timevar)) {
      df_grid <- data.frame(rep(data[[j]], each = length(grid[[1]])))
      df <- data.frame(data[[j]])
      names(df_grid) <- names(df) <- j
      X <- if(is.null(X)) df_grid else cbind(X, df_grid)
      XT <- if(is.null(XT)) df else cbind(XT, df)
    }
  }
  if(!is.null(X)) {
    colnames(X) <- colnames(XT) <- x$term[!(x$term %in% c(yname, timevar))]
  }

  X <- if(is.null(X)) data.frame(unlist(grid)) else cbind(X, unlist(grid))
  XT <- if(is.null(XT)) data.frame(y[, 1]) else cbind(XT, y[, 1])
  colnames(X)[ncol(X)] <- colnames(XT)[ncol(XT)] <- yname
  if(timevar != yname) {
    X <- cbind(X, unlist(grid))
    XT <- cbind(XT, y[, 1])
    colnames(X)[ncol(X)] <- colnames(XT)[ncol(XT)] <- timevar
  }
  if(x$by != "NA" & x$by != yname) {
    X[[x$by]] <- rep(data[[x$by]], each = length(grid[[1]]))
    XT[[x$by]] <- data[[x$by]]
    if (nrow(XT) != nrow(data)) {
      stop("XT dimensions do not coincide with 'by' dimensions for ", x$label)
    }
  }

  x$Xgrid <- mgcv::PredictMat(x, X)

  # XT necessary for calculation of score and hessian
  # Matrix of evaluations at the vector of survival times
  x$XT <- mgcv::PredictMat(x, XT)

  x
}



# PCRE Transformer Function -----------------------------------------------

sm_time_transform_mjm_pcre <- function(x, data, grid, yname, timevar, take,
                                       nmarker, cp_info) {
  if(!is.null(take))
    data <- data[take, , drop = FALSE]
  X <- NULL
  for(j in x$term) {
    if((j != yname) & (j != timevar)) {
      df <- data.frame(rep(data[[j]], each = length(grid[[1]])))
      names(df) <- j
      X <- if(is.null(X)) df else cbind(X, df)
    }
  }
  if(!is.null(X))
    colnames(X) <- x$term[!(x$term %in% c(yname, timevar))]
  X <- if(is.null(X)) data.frame(unlist(grid)) else cbind(X, unlist(grid))
  colnames(X)[ncol(X)] <- yname
  if(timevar != yname) {
    X <- cbind(X, unlist(grid))
    colnames(X)[ncol(X)] <- timevar
    data[[timevar]] <- data[[yname]]
  }
  if(x$by != "NA" & x$by != yname)
    X[[x$by]] <- rep(data[[x$by]], each = length(grid[[1]]))

  if(!"unc_pcre.random.effect" %in% class(x)) {
    class(x) <- c("pcre2.random.effect", "pcre.random.effect", "random.effect")
  }
  x$term <- c(x$term, timevar)
  x$timevar <- timevar
  x$Xgrid <- PredictMat(x, X, n = nmarker*nrow(X))
  x$XT <- PredictMat(x, data, n = nmarker*nrow(data))
  x$cp_info <- cp_info

  x
}


# Linear Design Transformer -----------------------------------------------

param_time_transform_mjm <- function(x, formula, data, grid, yname, timevar,
                                     take, idvar, y, timevar2 = NULL) {

  X <- Xn <- tvar <- NULL

  # For time-varying covariates in lambda predictor (idvar is not NULL)
  if (!is.null(idvar)) {
    id <- data[[idvar]]
    for (j in names(data)) {
      if ((!grepl("Surv(", j, fixed = TRUE) &
           !grepl("Surv2(", j, fixed = TRUE)) & (j != yname) & (j != timevar)) {

        # split data per subject
        idata <- split(data[[j]], id)
        # check if timevarying variable
        temp <- lapply(1:length(idata), function(i){
          length(unique(idata[[i]])) > 1
        })
        if(any(unlist(temp))){
          tvar <- c(tvar, j)
          # extract unique time-varying values
          values <- lapply(1:length(idata), function(i){unique(idata[[i]])})
          # extract break points
          breaks <- lapply(1:length(idata), function(i){
            split(data[[timevar2]], id)[[i]][c(TRUE, diff(idata[[i]]) != 0)]})
          # transfer break points to evaluation grid
          igrid <- lapply(1:length(idata), function(i){
            if(length(breaks[[i]]) > 1){
              g <- cut(grid[[i]], breaks[[i]], labels=FALSE,
                       include.lowest = TRUE)
              g[is.na(g)] <- max(g, na.rm=TRUE) + 1
              g
            } else {
              rep(1, length(grid[[i]]))
            }})
          # evaluate var iable on that grid
          evalgrid <- lapply(1:length(idata), function(i){
            values[[i]][igrid[[i]]]})
          df <- data.frame(unlist(evalgrid))
          names(df) <- j
          X <- if (is.null(X))
            df
          else cbind(X, df)
          Xn <- c(Xn, j)
        }
      }
    }
  }
  if(!is.null(take))
    data <- data[take, , drop = FALSE]

  for(j in names(data)) {
    if((!grepl("Surv(", j, fixed = TRUE) & !grepl("Surv2(", j, fixed = TRUE)) &
       (j != yname) & (j != timevar) & !(j %in% tvar)) {
      df <- data.frame(rep(data[[j]], each = length(grid[[1]])))
      names(df) <- j
      X <- if(is.null(X)) df else cbind(X, df)
      Xn <- c(Xn, j)
    }
  }
  if(!is.null(X)) {
    colnames(X) <- Xn
  }
  X <- if(is.null(X)) data.frame(unlist(grid)) else cbind(X, unlist(grid))
  colnames(X)[ncol(X)] <- yname
  if(timevar != yname) {
    X <- cbind(X, unlist(grid))
    colnames(X)[ncol(X)] <- timevar
  }

  x$Xgrid <- model.matrix(formula, data = X)
  # If the timevariate in the longitudinal part is included, then change
  # longitudinal time to survival time
  if (grepl(timevar, deparse(formula[[2]]))) {
    form <- gsub(timevar, yname, deparse(formula[[2]]))
    formula <- update.formula(formula, as.formula(paste("~", form)))
  }
  x$XT <- model.matrix(formula, data = data)

  x
}
