#' Fit IDR to training data
#'
#' @description Fits isotonic distributional regression (IDR) to a training
#'   dataset with respect to the stochastic order (SD) order.
#'
#' @usage idrsd(y, X, pars = osqpSettings(verbose = FALSE, eps_abs =
#'   1e-5, eps_rel = 1e-5, max_iter = 10000L), progress = TRUE)
#'
#' @param y numeric vector (the response variable).
#' @param X object that fits to specification in \emph{type}. For \emph{type = 'idr'}, X must be an IDR object,
#' for \emph{type = 'ensemble'}, X must be a data.frame, for \emph{type = 'ecdf'}, X must be a matrix,
#' for \emph{type = 'dis'} X is empty, for \emph{type = 'normal'} and \emph{type = 'normal_ab'} X is a matrix with 2 columns
#' which represent mu and sigma parameters of the normal distribution
#' @param type default is 'ensemble'. Other possibilities are 'idr', 'ecdf', 'dis', 'normal', normal_ab
#' @param grid if type == 'ecdf', than grid is vector of threshold values corresponding to ECDF-values in X
#' @param dis_func if type == 'dis', then a cumulative distribution function must be specified along with its distributional parameters.
#' @param pars parameters for quadratic programming optimization (only relevant
#'   if \code{X} has more than one column), set using
#'   \code{\link[osqp]{osqpSettings}}.
#' @param progress display progressbar (\code{TRUE}, \code{FALSE} or \code{1},
#'   \code{0})?
#'
#' @details This function computes the isotonic distributional regression (IDR)
#'   of a response \emph{y} on distributional data \emph{X}. IDR estimates
#'   the cumulative distribution function (CDF) of \emph{y} conditional on
#'   \emph{X} under the assumption of isotonicity. The odering on \emph{X}
#'   is computed using the stochastic order.
#'   The conditional CDFs are estimated at each threshold in \code{unique(y)}.
#'   This is the set where the CDFs may have jumps.
#'
#' @return An object of class \code{"idrcal"} containing the following
#'   components:
#'
#'   \item{\code{X}}{data frame of all distinct values (note the different order!).}
#'
#'   \item{\code{y}}{list of all observed responses in the training data for
#'   given \code{X}.}
#'
#'   \item{\code{cdf}}{matrix containing the estimated CDFs, one CDF per row,
#'   evaluated at \code{thresholds} (see next point). The CDF in the \code{i}th
#'   row corredponds to the estimated conditional distribution of the response
#'   given \code{X[i,]}.}
#'
#'   \item{\code{thresholds}}{the thresholds at which the CDFs in \code{cdf} are
#'   evaluated. The entries in \code{cdf[,j]} are the conditional CDFs evaluated
#'   at \code{thresholds[j]}.}
#'
#'   \item{\code{type}}{type of distributional input data as specified by user}
#'
#'   \item{\code{diagnostic}}{list giving a bound on the precision of the CDF
#'   estimation (the maximal downwards-step in the CDF that has been detected)
#'   and the fraction of CDF estimations that were stopped at the iteration
#'   limit \code{max_iter}. Decrease the parameters \code{eps_abs} and/or
#'   \code{eps_rel} or increase \code{max_iter} in \code{pars} to improve the
#'   precision. See \code{\link[osqp]{osqpSettings}} for more optimization
#'   parameters.}
#'
#'   \item{\code{indices}}{ the row indices of  \code{X} in the
#'   original training dataset (used for in-sample predictions with
#'   \code{\link{predict.idrfit}}).}
#'
#'   \item{\code{constraints}}{ (in multivariate IDR, \code{NULL} otherwise)
#'   matrices giving the order constraints for optimization. Used in
#'   \code{\link{predict.idrfit}}.}
#'
#'
#' @note The function \code{idrsd} is only intended for fitting IDR model for a
#'   training dataset and storing the results for further processing, but not
#'   for prediction or evaluation, which is done using the output of
#'   \code{\link{predict.idrcal}}.
#'
#' @seealso The S3 method \code{\link{predict.idrcal}} for predictions based on
#'   an IDR fit.
#'
#' @export
#' @importFrom osqp osqp
#' @importFrom stats setNames
#' @importFrom scoringRules crps_sample
#' @importFrom scoringRules crps_norm
#' @importFrom scoringRules crps_cnorm


idrsd <- function(y, X=NULL, grid = NULL, dis_func = NULL, type = 'ensemble', # 'dis', 'ecdf', 'normal', 'idr', normal_ab
                  inta = NULL, intb = NULL, org_crps = FALSE,
                  pars = osqpSettings(verbose = FALSE, eps_abs = 1e-5,
                                     eps_rel = 1e-5, max_iter = 10000L),
                  progress = TRUE, ...) {

  # Check input
  if (!is.vector(y, mode = "numeric"))
    stop("'y' must be a numeric vector")
  thresholds <- sort(unique(y))
  if ((nThr <- length(thresholds)) == 1)
    stop("'y' must contain more than 1 distinct value")
  if (isTRUE(progress == 1)) progress <- TRUE
  if (isTRUE(progress == 0)) progress <- FALSE
  if (!isTRUE(progress) & !isFALSE(progress))
    stop("'progress' must be TRUE/FALSE or 1/0")


  if (isTRUE(progress == 1)) progress <- TRUE
  if (isTRUE(progress == 0)) progress <- FALSE
  if (!isTRUE(progress) & !isFALSE(progress))
    stop("'progress' must be TRUE/FALSE or 1/0")

  original_crps = NULL

  if (type == 'idr'){
    #print('check if cdf = 0 needs to be included?')
    if (class(X) == 'idr'){
      if (org_crps) {
        original_crps <- mean(crps(X, y))
      }
      grid <- lapply(X, function(x) x$points)
      #grid_unique <- sort(unique(unlist(threshold_list, use.names=FALSE)))
      X <- lapply(X, function(x) x$cdf)
      #type = 'ecdf'
    } else
      stop("wrong 'X' for type 'idr'")
  }

  if (type == 'ensemble'){
    if (!is.data.frame(X))
      stop("'X' must be a data.frame")
    if (!all(sapply(X, function(col) is.numeric(col) || is.ordered(col))))
      stop("'X' must contain numeric or ordered factor variables")
    if (nrow(X) <= 1)
      stop("'X' must have more than 1 rows")
    if (nrow(X) != length(y))
      stop("length(y) and nrow(X) must match")
    if (anyNA(X) || anyNA(y))
      stop("'X' and 'y' must not contain NAs")

    nVar <- ncol(X)
    oldNames <- names(X)
    names(X) <- NULL

    if (org_crps){
      original_crps <- mean(crps_sample(y, dat = X))
    }

    x <- data.frame(y = y, ind = seq_along(y))
    X <- data.frame(t(apply(X, 1, FUN=function(x) sort(x)))) # X can be sorted rowwise, since
    # the ensemble members are considered to be exchangeable.
    X <- stats::aggregate(x = x, by = X, FUN = identity, simplify = FALSE)
    cpY <- X[["y"]]
    indices <- X[["ind"]]
    X <- X[, 1:nVar, drop = FALSE] # compact X: no row twice
    names(X) <- oldNames
    weights <- sapply(indices, length)

    constr <- compOrd_cpp(X)
    N <- nrow(X)

  } else if (type == 'normal') {
    if (ncol(X) != 2)
      stop("'X' must contain two columns with mean and std")
    #if (!is.data.frame(X))
    # stop("'X' must be a data.frame")
    if (!all(sapply(X, function(col) is.numeric(col))))
      stop("'X' must contain numeric variables")
    if (anyNA(X) || anyNA(y))
      stop("'X' and 'y' must not contain NAs")
    if (nrow(X) != length(y))
      stop("length(y) and nrow(X) must match")
    if (!all(X[, 2]>=0))
      stop("sigma bust be positive")
    #if (eps == 'sd')
    #  stop("'eps = sd' not yet implemented for type = normal")
    if (org_crps) {
      original_crps <- mean(crps_norm(y, X[,1], X[,2]))
    }
    x <- data.frame(y = y, ind = seq_along(y))
    X <- stats::aggregate(x = x, by = data.frame(X), FUN = identity, simplify = FALSE)
    cpY <- X[["y"]]
    indices <- X[["ind"]]
    X <- data.matrix(X[, c(1, 2)], rownames.force = NA)
    #print(X[1:10,])
    constr <- compOrd_normal(X)
    weights <- sapply(indices, length)
    N <- nrow(X)
    nVar <- 2
  }  else if (type == 'normal_ab') {
    if (ncol(X) != 2)
      stop("'X' must contain two columns with mean and std")
    #if (!is.data.frame(X))
    # stop("'X' must be a data.frame")
    if (!all(sapply(X, function(col) is.numeric(col))))
      stop("'X' must contain numeric variables")
    if (anyNA(X) || anyNA(y))
      stop("'X' and 'y' must not contain NAs")
    if (nrow(X) != length(y))
      stop("length(y) and nrow(X) must match")
    if (!all(X[, 2]>=0))
      stop("sigma bust be positive")
    if (is.null(inta) || is.null(intb)){
      stop("define a and b")
    }
    if (inta >= intb){
      stop("a must be smaller than b")
    }

    if (org_crps) {
      original_crps <- mean(crps_cnorm(y, X[,1], X[,2], lower = inta, upper = intb))
    }
    #if (eps == 'sd')
    #  stop("'eps = sd' not yet implemented for type = normal")
    x <- data.frame(y = y, ind = seq_along(y))
    X <- stats::aggregate(x = x, by = data.frame(X), FUN = identity, simplify = FALSE)
    cpY <- X[["y"]]
    indices <- X[["ind"]]
    X <- data.matrix(X[, c(1, 2)], rownames.force = NA)
    #print(X[1:10,])
    constr <- compOrd_normal_ab(X, inta, intb)
    weights <- sapply(indices, length)
    N <- nrow(X)
    nVar <- 2
  } else if (type == 'dis') {
    z <- list(...)
    para_dim <- length(z)
    #if ( eps == 'sd')
     # stop("'eps = sd' not yet implemented")
    if (para_dim == 0)
      stop("arguments for dis_func must be non empty")
    #if (!all(sapply(z, function(x) length(x) == length(y))))
     # stop("length of each argument in dis_func and length(y) must match")
    if (is.null(dis_func))
      stop("'dis_func' must be non-empty")
    if (is.null(grid))
      stop("'grid' must be non-empty")
    if (!is.vector(grid, mode = "numeric"))
      stop("'grid' must be a numeric vector")

    nVar <- length(grid)
    list_names = names(z)
    P <- matrix(unlist(z), ncol =para_dim , nrow = length(y))
    P <- as.data.frame(P)
    colnames(P) <- list_names
    x <- data.frame(y = y, ind = seq_along(y))
    P2 <- stats::aggregate(x = x, by = P, FUN = identity, simplify = FALSE)
    cpY <- P2[["y"]]
    indices <- P2[["ind"]]
    #print(P2[1:10, c(1,2)])
    #P2 <- as.matrix(P2[, seq(1, para_dim, 1)], rownames.force = NA)
    #P <- as.list(P2[, seq(1, para_dim, 1)])
    d0 <- length(grid)
    #first <- dis_func(grid[1], ...)
    n0 <- length(y)
    X <- matrix(0, nrow = n0, ncol= d0)

    for(i in 1:d0){
      X[,i] <- dis_func(grid[i], ...)
    }

    if (org_crps) {
      original_crps <- mean(ecdf_crps(y, grid, X))
    }

    # ECDF in X with grid
    indices_sel <-  sapply(indices, FUN = function(x) x[1])
    X <- X[indices_sel, ]

    constr <- compOrd_ecdf(X, grid)
    weights <- sapply(indices, length)
    N <- length(cpY)

  } else if (type == 'ecdf' || type == 'idr') {
    if (is.vector(grid, mode = "numeric") && (is.matrix(X))){
      if (length(grid) != ncol(X))
        stop("length(grid) and ncol(X) must match")
      if (nrow(X) != length(y))
        stop("length(y) and nrow(X) must match")
      if (anyNA(X) || anyNA(y))
        stop("'X' and 'y' must not contain NAs")
      if (anyNA(grid))
        stop("'grid' must not contain NAs")
    } else if (is.list(grid) && is.list(X)) {
      if (length(X) != length(y))
        stop('length(y) and length(X) must match')
      if (length(grid) != length(y))
        stop('length(y) and length(grid) must match')
      length_X <- sapply(X, function(x) length(x))
      length_grid <- sapply(grid, function(x) length(x))
      if (!all(length_X == length_grid))
        stop("length of list elements in X and grid must match")
    } else {
      stop("invalid values of 'X' and 'grid'")
    }

    if (org_crps) {
      if (type == 'ecdf'){
        original_crps <- mean(ecdf_crps(y, grid, X))
      }
    }
    if (is.list(grid)){
      #M_class <- ecdf_list_comp_class_sd(X, grid)
      stop("'grid' must be a vector")
    } else {
      M_class <- ecdf_comp_class_sd(X, grid)
    }

    class_X <- M_class$ind


    M <- M_class$M

    x <- data.frame(y = y, ind = seq_along(y))

    split_classes <- stats::aggregate(x = x, by = data.frame(class_X), FUN = identity, simplify = FALSE)

    cpY <- split_classes[["y"]]
    indices <- split_classes[["ind"]]
    class_X_sorted <- split_classes[, 1]
    vec_indx <- sapply(indices, FUN = function(x) x[1])#unlist(indices)
    M <- M[vec_indx, vec_indx]
    paths <- which(M==1,arr.ind = TRUE)
    dimnames(paths) <- list(NULL, c("smaller", "greater"))
    constr <- list(M=M, paths = paths)
    weights <- sapply(indices, length)
    nVar <- length(class_X)
    N <- length(cpY)

    if (is.list(X)){
      X <- sapply(vec_indx, function(i) X[[i]] )
      grid <- sapply(vec_indx, function(i) grid[[i]] )
    } else {
      X <- X[vec_indx, ]
    }
    #print(weights)
    #X <- data.matrix(X2[, c(1, 2)], rownames.force = NA)
    #constr <- normal_comp(X, eps)
    #weights <- sapply(indices, length)
    #N <- nrow(X)

    # List X and List grid

    # Two options ?

  } else {
    stop("invalid value for 'type'")
  }

  # should work with normal and ensemble
  cdf <- matrix(ncol = nThr - 1, nrow = N)
  A <- trReduc(constr$paths, N)
  nConstr <- nrow(A)
  l <- rep(0, nConstr)
  A <- Matrix::sparseMatrix(i = rep(seq_len(nConstr), 2), j = as.vector(A),
                            x = rep(c(1, -1), each = nConstr), dims = c(nConstr, N))


  P <- Matrix::sparseMatrix(i = 1:N, j = 1:N, x = weights)
  i <- 1
  I <- nThr - 1
  conv <- vector("logical", I)
  q <- - weights * sapply(cpY, FUN = function(x) mean(thresholds[i] >= x))
  qp <- osqp::osqp(P = P, q = q, A = A, l = l, pars = pars)
  sol <- qp$Solve()
  cdf[, 1] <- pmin(1, pmax(0, sol$x))
  conv[1] <- identical(sol$info$status, "maximum iterations reached")

  if (I > 1) {
    if (progress) {
      cat("Estimating cdf...\n")
      pb <- utils::txtProgressBar(style = 1)
      for (i in 2:I) {
        utils::setTxtProgressBar(pb, i/I)
        qp$WarmStart(x = cdf[, i - 1L])
        q <-  -weights * sapply(cpY, FUN = function(x) mean(thresholds[i] >= x))
        qp$Update(q = q)
        sol <- qp$Solve()
        cdf[, i] <- pmin(1, pmax(0, sol$x))
        conv[i] <- identical(sol$info$status, "maximum iterations reached")
      }
      close(pb)
      cat("\n")
    } else {
      for (i in 2:I) {
        qp$WarmStart(x = cdf[, i - 1L])
        q <-  -weights * sapply(cpY, FUN = function(x) mean(thresholds[i] >= x))
        qp$Update(q = q)
        sol <- qp$Solve()
        cdf[, i] <- pmin(1, pmax(0, sol$x))
        conv[i] <- identical(sol$info$status, "maximum iterations reached")
      }
    }
  }
  diagnostic <- list(
    precision = ifelse(I > 1, abs(min(diff(t(cdf)))), 0),
    convergence = mean(conv)
  )

  # Apply pava to estimated CDF to ensure monotonicity
  if (nVar > 1) cdf <- cbind(pavaCorrect(cdf), 1)


  # How to define output ?
  if (type == 'ecdf' || type == 'dis' || type == 'idr'){
    structure(list(X = X, y = cpY, cdf = cdf, thresholds = thresholds, type = type, grid = grid,
                   diagnostic = diagnostic, indices = indices, Ord=constr, org_crps = original_crps),
              class = "idrcal")
  } else {
    structure(list(X = X, y = cpY, cdf = cdf, thresholds = thresholds, type = type,
                   diagnostic = diagnostic, indices = indices, Ord=constr, org_crps = original_crps),
              class = "idrcal")
  }



}


#' Predict method for IDR fits
#'
#' @description Prediction based on IDR model fit.
#'
#' @method predict idrcal
#'
#' @param object IDR fit (object of class \code{"idrcal"}).
#' @param data optional object with which to
#'   predict. In-sample predictions are returned if this is omitted.
#' @param digits number of decimal places for the predictive CDF.
#' @param ... included for generic function consistency.
#'
#' @details If the variables \code{data} for which predictions are
#' desired are already contained in the training dataset \code{X} for the fit,
#' \code{predict.idrfit} returns the corresponding in-sample prediction.
#' Otherwise monotonicity is used to derive upper and lower bounds for the
#' predictive CDF, and the predictive CDF is a pointwise average of these
#' bounds.
#'
#' @return A list of predictions. Each prediction is a \code{data.frame}
#'   containing the following variables:
#'
#' \item{\code{points}}{the points where the predictive CDF has jumps.}
#'
#' \item{\code{cdf}}{the estimated CDF evaluated at the \code{points}.}
#'
#' \item{\code{lower}, \code{upper}}{ (only for out-of-sample predictions)
#'   bounds for the estimated CDF, see 'Details' above.}
#'
#' The output has the attribute \code{incomparables}, which gives the indices
#' of all predictions for which the climatological forecast is returned because
#' the forecast variables are not comparable to the training data.
#'
#' @export
#' @importFrom stats predict
#' @importFrom stats approx
#'
#' @seealso
#' \code{\link{idrsd}} to fit IDR to training data.
#'
#' \code{\link{cdf}}, \code{\link{qpred}} to evaluate the CDF or quantile
#' function of IDR predictions.
#'
#' \code{\link{bscore}}, \code{\link{qscore}}, \code{\link{crps}},
#' \code{\link{pit}} to compute Brier scores, quantile scores, the CRPS and the
#' PIT of IDR predictions.
#'
#' \code{\link[isodistrreg:plot.idrsd]{plot}} to plot IDR predictive CDFs.
#'
#'


predict.idrcal <- function(object, data = NULL, grid = NULL, digits = 3,
                           interpolation = "linear", ...) {
  cdf <- object$cdf
  thresholds <- object$thresholds

  z <- list(...)

  if (is.null(data) && length(z) == 0) {
    # In-sample predictions
    indices <- object$indices
    preds <- structure(vector("list", length(unlist(indices))),
                       class = c("idrsd"))
    for (i in seq_along(indices)) {
      edf <- round(cdf[i, ], digits)
      sel <- c(edf[1] > 0, diff(edf) > 0)
      tmp <- data.frame(points = thresholds[sel], cdf = edf[sel])
      for (j in indices[[i]]) {
        preds[[j]] <- tmp
      }
    }
    return(preds)
  }

  # Out-of-sample predictions
  type <- object$type
  X <- object$X

  if (type == 'idr'){
    if (class(data) == 'idr'){
      # if matrix:
      #points <- lapply(data, function(x) x$points)
      #grid <- sort(unique(unlist(points, recursive = TRUE, use.names = FALSE)))
      #data <- lapply(data, function(x) x$cdf)
      #data <- ecdf_func(data, points,  grid)

      # if list:
      grid <- lapply(data, function(x) x$points)
      data <- lapply(data, function(x) x$cdf)
      #type = 'ecdf'
    } else
      stop("wrong 'data' for type 'idr'")
  }
    # prepare validation data and run ecdf


  if (type == 'ecdf' || type == 'idr'){
    grid_base <- object$grid

    nPoints <- neighborPoints3(x = data, X = X, gridx = grid, gridX = grid_base, Ord = object$Ord)

    if (is.matrix(data)){
      nx = nrow(data)
    } else {
      nx = length(data)
    }
  }

  else if (type == 'ensemble'){
    if (!is.data.frame(data))
      stop("'data' must be a data.frame")
    M <- match(colnames(X), colnames(data), nomatch = 0)
    if (any(M == 0))
      stop("some variables of the idr fit are missing in 'data'")
    nVar <- ncol(data)
    nx <- nrow(data)

    # Prediction Method for multivariate IDR
    nPoints <- neighborPoints(x = data, X = X, Ord = object$Ord)
    # run predict()
  }

  else if (type == 'dis'){
    grid_base <- object$grid
    d0 <- length(grid_base)
    #first <- dis_func(grid[1], ...)
    n0 <- length(z[[1]])
    data <- matrix(0, nrow = n0, ncol= d0)

    for(i in 1:d0){
      data[,i] <- dis_func(grid_base[i], ...)
    }

    # neighbourPoints smarter implementation since defined over same set of thresholds
    nPoints <- neighborPoints_samegrid(x = data, X = X, grid = grid_base, Ord = object$Ord)
    nx <- nrow(data)
  }

  else if (type == 'normal'){
    #if (eps == 'sd') {
    #  stop("'eps = sd' not yet implemented for type = normal")
    #}
    #else {
    nPoints <- neighborPoints_norm(x = data, X = X, Ord = object$Ord)
    nx <- nrow(data)
    #}

  }

  #if (!is.data.frame(data))
   # stop("'data' must be a data.frame")

  #M <- match(colnames(X), colnames(data), nomatch = 0)
  #if (any(M == 0))
  #  stop("some variables of the idr fit are missing in 'data'")
  #nVar <- ncol(data)
  #nx <- length(data)

  # Prediction Method for multivariate IDR

  preds <- structure(vector("list", nx), class = c("idrsd"))
  smaller <- nPoints$smaller
  greater <- nPoints$greater


  # Climatological forecast for incomparable variables
  incomparables <- sapply(smaller, length) + sapply(greater, length) == 0
  if (any(incomparables)) {
    y <- unlist(object$y)
    edf <- round(stats::ecdf(y)(thresholds), digits)
    sel <- edf > 0
    edf <- edf[sel]
    points <- thresholds[sel]
    upr <- which.max(edf == 1)
    if (upr < length(edf)) {
      points <- points[-((upr + 1):length(edf))]
      edf <- edf[-((upr + 1):length(edf))]
    }
    dat <- data.frame(points = points, lower = edf, cdf = edf, upper = edf)
    for (i in which(incomparables)) preds[[i]] <- dat
  }

  # Predictions for comparable variables
  for (i in which(!incomparables)) {
    # Case distinction: Existence of lower and/or upper bound on CDF
    if (length(smaller[[i]]) > 0 && length(greater[[i]]) == 0) {
      upper <- round(apply(cdf[smaller[[i]], , drop = FALSE], 2, min), digits)
      sel <- c(upper[1] != 0, diff(upper) != 0)
      upper <- estimCdf <- upper[sel]
      lower <- rep(0, length(upper))
    } else if (length(smaller[[i]]) == 0 && length(greater[[i]]) > 0) {
      lower <- round(apply(cdf[greater[[i]], , drop = FALSE], 2, max), digits)
      sel <- c(lower[1] != 0, diff(lower) != 0)
      lower <- estimCdf <- lower[sel]
      upper <- rep(1, length(lower))
    } else {
      lower <- round(apply(cdf[greater[[i]], , drop = FALSE], 2, max), digits)
      upper <- round(apply(cdf[smaller[[i]], , drop = FALSE], 2, min), digits)

      sel <- c(lower[1] != 0, diff(lower) != 0) |
        c(upper[1] != 0, diff(upper) != 0)

      lower <- lower[sel]
      upper <- upper[sel]
      estimCdf <- round((lower + upper) / 2, digits)
    }
    preds[[i]] <- data.frame(points = thresholds[sel], lower = lower,
                             cdf = estimCdf, upper = upper)
  }
  attr(preds, "incomparables") <- which(incomparables)
  preds
}


#' @export
print.idrsd <- function(x, ...) {
  # Print only head of first prediction
  cat(paste("IDR predictions: list of", length(x), "prediction(s)\n\n"))
  print(list(utils::head(x[[1]])))
  cat("...")
  invisible(x)
}


#' @export
print.idrcal <- function(x, ...) {
  # Print diagnostic information
  cat("IDR fit: \n")
  cat(paste("CDFs estimated:", nrow(x$cdf), "\n"))
  cat(paste("Thresholds for estimation:", ncol(x$cdf), "\n"))
  prec <- signif(x$diagnostic$precision, 2)
  conv <- signif(x$diagnostic$convergence, 4) * 100
  cat(paste("CDF estimation precision:", prec, "\n"))
  cat(paste0("Estimations stopped after max_iter iterations: ", conv, "% \n"))
  invisible(x)
}

#' @importFrom osqp osqpSettings
#' @export
osqp::osqpSettings
