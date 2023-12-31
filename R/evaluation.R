#' Cumulative distribution function (CDF) of IDRsd or raw forecasts
#'
#' @description
#' Evaluate the the cumulative distribution function (CDF) of IDRsd predictions or
#' of unprocessed forecasts in a \code{data.frame}.
#'
#' @usage
#' cdf(predictions, thresholds)
#'
#' @param predictions either an object of class \code{idrsd} (output of
#'   \code{\link{predict.idrcal}}), or a \code{data.frame} of numeric variables.
#'   In the latter case, the CDF is computed using the empirical distribution of
#'   the variables in \code{predictions}.
#' @param thresholds numeric vector of thresholds at which the CDF will be
#'   evaluated.
#'
#' @details
#' The CDFs are considered as piecewise constant stepfunctions: If \code{x} are
#' the points where the IDR fitted CDF (or the empirical distribution of the
#' forecasts) has jumps and \code{p} the corresponding CDF values, then for
#' \code{x[i] <= x < x[i + 1]}, the CDF at \code{x} is \code{p[i]}.
#'
#' @return
#' A matrix of probabilities giving the evaluated CDFs at the given thresholds,
#' one column for each threshold.
#'
#' @seealso
#' \code{\link{predict.idrcal}} \code{\link{qpred}}, \code{\link{bscore}}
#'
#' @export
#'
#' @examples
#' ## Data from IDR package:
#' library(isodistrreg)
#' data("rain")
#'
#' ## IDRsd based on ensemble forecast
#'
#' ensemble <- rain[1:(3 * 365), 3:54, drop = FALSE]
#' y <- rain[1:(3 * 365), "obs"]
#'
#' fit <- idrsd(y = y, X = ensemble, type = 'ensemble')
#' predictions <- predict(fit)
#'
#' ## Compute probability of precipitation
#' 1 - cdf(predictions, thresholds = 0)
#'
cdf <- function(predictions, thresholds) {
    UseMethod("cdf")
}

#' cdf method for class 'idr'
#'
#' @method cdf idrsd
#' @rdname cdf
#' @importFrom stats stepfun
#' @export
cdf.idrsd <- function(predictions, thresholds) {
  if (!is.vector(thresholds, "numeric"))
    stop("'thresholds' must be a numeric vector")
  cdf0 <- function(data) {
    # Evaluate CDF (stepfun) at thresholds
    stats::stepfun(x = data$points, y = c(0, data$cdf))(thresholds)
  }
  cdfVals <- lapply(predictions, cdf0)
  do.call(rbind, cdfVals)
}

#' cdf method for class 'data.frame'
#'
#' @method cdf data.frame
#' @rdname cdf
#' @importFrom stats ecdf
#' @export
cdf.data.frame <- function(predictions, thresholds) {
  if (!is.vector(thresholds, "numeric"))
    stop("'thresholds' must be a numeric vector")
  if (!all(sapply(predictions, is.numeric)))
      stop("'predictions' contains non-numeric variables")
  n <- nrow(predictions)
  predictions <- unname(split(data.matrix(predictions), seq_len(n)))
    # Update to asplit instead of split(data.matrix(...))
  cdf0 <- function(data) stats::ecdf(data)(thresholds)
  cdfVals <- lapply(predictions, cdf0)
  do.call(rbind, cdfVals)
}


#' Quantile function of IDRsd or raw forecasts
#'
#' @description
#' Evaluate the the quantile function of IDRsd predictions or of unprocessed
#' forecasts in a \code{data.frame}.
#'
#' @usage
#' qpred(predictions, quantiles)
#'
#' @param predictions either an object of class \code{idrsd} (output of
#'   \code{\link{predict.idrcal}}), or a \code{data.frame} of numeric variables. In
#'   the latter case, quantiles are computed using the empirical distribution of
#'   the variables in \code{predictions}.
#' @param quantiles numeric vector of desired quantiles.
#'
#' @details
#' The quantiles are defined as lower quantiles, that is,
#' \deqn{
#'   q(u) = inf(x: cdf(x) >= u).
#' }
#'
#' @return
#' A matrix of forecasts for the desired quantiles, one column per quantile.
#'
#' @seealso
#' \code{\link{predict.idrcal}}, \code{\link{cdf}}, \code{\link{qscore}}
#'
#' @export
#'
#' @examples
#' ## Data from IDR package:
#' library(isodistrreg)
#' data("rain")
#'
#' ## IDRsd based on ensemble forecast
#'
#' ensemble <- rain[1:(3 * 365), 3:54, drop = FALSE]
#' y <- rain[1:(3 * 365), "obs"]
#'
#' fit <- idrsd(y = y, X = ensemble, type = 'ensemble')
#' predictions <- predict(fit)
#'
#' ## Compute 95%-quantile forecast
#' qpred(predictions, quantiles = 0.95)
#'
qpred <- function(predictions, quantiles) {
    UseMethod("qpred")
}

#' qpred method for class 'idrsd'
#'
#' @method qpred idrsd
#'
#' @importFrom stats stepfun
#' @importFrom utils tail
#' @rdname qpred
#' @export
qpred.idrsd <- function(predictions, quantiles) {

  # Check input
  if (!is.vector(quantiles, "numeric") || min(quantiles) < 0 ||
      max(quantiles) > 1)
    stop("quantiles must be a numeric vector with entries in [0,1]")
  q0 <- function(data) {
    # Evaluate quantile function (stepfun) at given quantiles
    stats::stepfun(x = data$cdf,
      y = c(data$points, data$points[nrow(data)]), right = TRUE)(quantiles)
  }
  qVals <- lapply(predictions, q0)
  do.call(rbind, qVals)
}


#' qpred method for class 'data.frame'
#'
#' @method qpred data.frame
#'
#' @rdname qpred
#' @importFrom stats quantile
#' @export
qpred.data.frame <- function(predictions, quantiles) {

  # Check input
  if (!is.vector(quantiles, "numeric") || min(quantiles) < 0 ||
      max(quantiles) > 1)
    stop("quantiles must be a numeric vector with entries in [0,1]")
  if (!all(sapply(predictions, is.numeric)))
      stop("'predictions' contains non-numeric variables")

  n <- nrow(predictions)
  predictions <- split(data.matrix(predictions), seq_len(n))
    # Update to asplit instead of split(data.matrix(...))
  q0 <- function(x) stats::quantile(x, probs = quantiles, type = 1)
  qVals <- lapply(predictions, q0)
  unname(do.call(rbind, qVals))
}

#' Quantile scores for IDRsd or raw forecasts
#'
#' @description
#' Computes quantile scores of IDRsd quantile predictions or of quantile
#' predictions from raw forecasts in a \code{data.frame}.
#'
#' @usage
#' qscore(predictions, quantiles, y)
#'
#' @inheritParams qpred
#' @param y a numeric vector of obervations of the same length as the number of
#'   predictions, or of length 1. In the latter case, \code{y} will be used
#'   for all predictions.
#'
#' @details
#' The quantile score of a forecast \emph{x} for the \emph{u}-quantile is
#' defined as
#' \deqn{
#' 2(1{x > y} - u)(x - y),
#' }
#' where \emph{y} is the observation. For \emph{u = 1/2}, this equals the mean
#' absolute error of the median forecast.
#'
#' @return
#' A matrix of the quantile scores for the desired quantiles, one column per
#' quantile.
#'
#' @seealso
#' \code{\link{predict.idrcal}}, \code{\link{qpred}}
#'
#' @references
#' Gneiting, T. and Raftery, A. E. (2007), 'Strictly proper scoring rules,
#' prediction, and estimation', Journal of the American Statistical Association
#' 102(477), 359-378
#'
#' @export
#'
#' @examples
#' data("rain")
#'
#' ## IDRsd based on ensemble forecast
#'
#' ensemble <- rain[1:(3 * 365), 3:54, drop = FALSE]
#' y <- rain[1:(3 * 365), "obs"]
#'
#' fit <- idrsd(y = y, X = ensemble, type = 'ensemble')
#'
#' ## Compute mean absolute error of the median postprocessed forecast using
#' ## data of the next year (out-of-sample predictions) and compare to raw
#' ## forecast
#'
#' ytest = rain[(3 * 365 + 1):(4 * 365), "obs"]
#' ensemble_test = rain[(3 * 365 + 1):(4 * 365), 3:54,drop = FALSE]
#' predictions <- predict(fit, data = ensemble_test)
#'
#' idrMAE <- mean(qscore(predictions, 0.5, ytest))
#' rawMAE <- mean(qscore(ensemble_test, 0.5, ytest))
#'
#' c("idrsd" = idrMAE, "raw" = rawMAE)
#'
qscore <- function(predictions, quantiles, y) {
  if (!is.vector(y, "numeric"))
    stop("y must be a numeric vector")
  predicted <- qpred(predictions, quantiles)
  if (!(lY <- length(y)) %in% c(1, nrow(predicted)))
    stop("y must have length 1 or the same length as the predictions")
  qsVals <- sweep(x = predicted, STATS = y, MARGIN = 1)
  2 * qsVals * sweep(qsVals > 0, MARGIN = 2, STATS = quantiles)
}

#' Brier score for forecast probability of threshold exceedance
#'
#' @description Computes the Brier score of forecast probabilities for exceeding
#' given thresholds.
#'
#' @usage bscore(predictions, thresholds, y)
#'
#' @inheritParams cdf
#' @param y a numeric vector of obervations of the same length as the number of
#'   predictions, or of length 1. In the latter case, \code{y} will be used for
#'   all predictions.
#'
#' @details The Brier score for the event of exceeding a given threshold
#' \emph{z} is defined as \deqn{ (1\{y > z\} - P(y > z))^2 } where \emph{y} is the
#' observation and \emph{P(y > z)} the forecast probability for exceeding the
#' threshold \code{z}.
#'
#' @return A matrix of the Brier scores for the desired thresholds, one column
#' per threshold.
#'
#' @seealso \code{\link{predict.idrcal}}, \code{\link{cdf}}
#'
#' @references
#' Gneiting, T. and Raftery, A. E. (2007), 'Strictly proper scoring rules,
#' prediction, and estimation', Journal of the American Statistical Association
#' 102(477), 359-378
#'
#' @export
#'
#' @examples
#' data("rain")
#'
#' ## IDRsd based on ensemble forecast
#'
#' ensemble <- rain[1:(3 * 365), 3:54, drop = FALSE]
#' y <- rain[1:(3 * 365), "obs"]
#'
#' fit <- idrsd(y = y, X = ensemble, type = 'ensemble')
#'
#' ## Compute Brier score for probability of precipitation
#' ## forecast using data of the next year (out-of-sample predictions)
#'
#' ytest = rain[(3 * 365 + 1):(4 * 365), "obs"]
#' ensemble_test = rain[(3 * 365 + 1):(4 * 365), 3:54,drop = FALSE]
#' predictions <- predict(fit, data = ensemble_test)
#' score <- bscore(predictions, thresholds = 0, y = ytest)
#'
#' mean(score)
bscore <- function(predictions, thresholds, y) {

  # Check input
  if (!is.vector(y, "numeric"))
    stop("obs must be a numeric vector")
  predicted <- cdf(predictions, thresholds)
  if (!(lY <- length(y)) %in% c(1, nrow(predicted)))
    stop("y must have length 1 or the same length as the predictions")

  if (lY == 1) {
    bsVals <- sweep(predicted, MARGIN = 2, STATS = (y <= thresholds))
  } else {
    bsVals <- predicted - outer(y, thresholds, "<=")
  }
  bsVals^2
}

#' Continuous ranked probability score (CRPS)
#'
#' @description Computes the CRPS of IDRsd or raw forecasts.
#'
#' @usage crps(predictions, y)
#'
#' @param predictions either an object of class \code{idrsd} (output of
#'   \code{\link{predict.idrcal}}), or a \code{data.frame} of numeric variables. In
#'   the latter case, the CRPS is computed using the empirical distribution of
#'   the variables in \code{predictions}.
#' @param y a numeric vector of obervations of the same length as the number of
#'   predictions, or of length 1. In the latter case, \code{y} will be used for
#'   all predictions.
#'
#' @details
#' This function uses adapted code taken from the function \code{crps_edf} of
#' the \pkg{scoringRules} package.
#'
#' @return A vector of CRPS values.
#'
#' @seealso \code{\link{predict.idrcal}}
#'
#' @references
#' Jordan A., Krueger F., Lerch S. (2018). "Evaluating Probabilistic
#' Forecasts with scoringRules." Journal of Statistical Software. Forthcoming.
#'
#' Gneiting, T. and Raftery, A. E. (2007), 'Strictly proper scoring rules,
#' prediction, and estimation', Journal of the American Statistical Association
#' 102(477), 359-378
#'
#' @export
#'
#' @examples
#' data("rain")
#'
#' ## IDRsd based on ensemble forecast
#'
#' ensemble <- rain[1:(3 * 365), 3:54, drop = FALSE]
#' y <- rain[1:(3 * 365), "obs"]
#'
#' fit <- idrsd(y = y, X = ensemble, type = 'ensemble')
#'
#' ## Compute CRPS of forecast using data of the next year
#' ## (out-of-sample predictions)
#'
#' ytest = rain[(3 * 365 + 1):(4 * 365), "obs"]
#' ensemble_test = rain[(3 * 365 + 1):(4 * 365), 3:54,drop = FALSE]
#' predictions <- predict(fit, data = ensemble_test)
#' idrCrps <- crps(predictions, y = ytest)
#'
#' ## Compare this to CRPS of the raw ensemble of all forecasts (high resolution,
#' ## control and 50 perturbed ensemble forecasts)
#'
#' rawData <- rain[(3 * 365 + 1):(4 * 365), c("HRES", "CTR", paste0("P", 1:50))]
#' rawCrps <- crps(rawData, y = ytest)
#'
#' c("idrsd" = mean(idrCrps), "raw_all" = mean(rawCrps))
crps <- function(predictions, y) {
    UseMethod("crps")
}

#' crps method for class 'idrsd'
#'
#' @method crps idrsd
#' @rdname crps
#' @export
crps.idrsd <- function(predictions, y) {

  # Check input
  if (!is.vector(y, "numeric"))
    stop("obs must be a numeric vector")
  if (length(y) != 1 && length(y) != length(predictions))
    stop("y must have length 1 or the same length as the predictions")

  x <- lapply(predictions, function(p) p$points)
  p <- lapply(predictions, function(p) p$cdf)
  w <- lapply(p, function(x) c(x[1], diff(x)))

  crps0 <- function(y, p, w, x) 2 * sum(w * ((y < x) - p + 0.5 * w) * (x - y))
  mapply(crps0, y = y, p, w = w, x = x)
}



#' crps method for class 'irdpred'
#'
#' @method crps idr
#' @rdname crps
#' @export
crps.idr <- function(predictions, y) {

  # Check input
  if (!is.vector(y, "numeric"))
    stop("obs must be a numeric vector")
  if (length(y) != 1 && length(y) != length(predictions))
    stop("y must have length 1 or the same length as the predictions")

  x <- lapply(predictions, function(p) p$points)
  p <- lapply(predictions, function(p) p$cdf)
  w <- lapply(p, function(x) c(x[1], diff(x)))

  crps0 <- function(y, p, w, x) 2 * sum(w * ((y < x) - p + 0.5 * w) * (x - y))
  mapply(crps0, y = y, p, w = w, x = x)
}



#' crps method for class 'data.frame'
#'
#' @method crps data.frame
#' @rdname crps
#' @export
crps.data.frame <- function(predictions, y) {
  # Check input
  if (!is.vector(y, "numeric"))
    stop("obs must be a numeric vector")
  n <- nrow(predictions)
  if (length(y) != 1 && length(y) != n)
    stop("y must have length 1 or the same length as the predictions")

  x <- lapply(split(data.matrix(predictions), seq_len(n)), sort)
  m <- ncol(predictions)
  a <- seq.int(0.5 / m, 1 - 0.5 / m, length.out = m)
  crps0 <- function(y, x, a) 2 / m * sum(((y < x) - a) * (x - y))
  mapply(crps0, y = y, x = x, MoreArgs = list(a = a))
}


#' Probability integral transform (PIT)
#'
#' @description Computes the probability integral transform (PIT) of IDRsd or raw
#'   forecasts.
#'
#' @usage pit(predictions, y, randomize = TRUE, seed = NULL)
#'
#' @param predictions either an object of class \code{idrsd} (output of
#'   \code{\link{predict.idrcal}}), or a \code{data.frame} of numeric variables. In
#'   the latter case, the PIT is computed using the empirical distribution of
#'   the variables in \code{predictions}.
#' @param y a numeric vector of obervations of the same length as the number of
#'   predictions.
#' @param randomize PIT values should be randomized at discontinuity points of the
#'   predictive CDF (e.g. at zero for precipitation forecasts). Set \code{
#'   randomize = TRUE} to randomize.
#' @param seed argument to \code{set.seed} for random number generation (if
#'   \code{randomize} is \code{TRUE}).
#'
#' @return Vector of PIT values.
#'
#' @seealso \code{\link{predict.idrcal}}
#'
#' @references
#'
#' Gneiting, T., Balabdaoui, F. and Raftery, A. E. (2007), 'Probabilistic
#' forecasts, calibration and sharpness', Journal of the Royal Statistical
#' Society: Series B (Statistical Methodology) 69(2), 243-268.
#'
#' @export
#'
#' @importFrom stats runif
#'
#' @examples
#' data("rain")
#' require("graphics")
#'
#'
#' ensemble <- rain[1:(3 * 365), 3:54, drop = FALSE]
#' y <- rain[1:(3 * 365), "obs"]
#'
#' fit <- idrsd(y = y, X = ensemble, type = 'ensemble')
#'
#' ## Assess forecast using data of next 2
#' ## years and compare to calibration of the raw ensemble
#'
#' ytest = rain[(3 * 365 + 1):(5 * 365), "obs"]
#' ensemble_test = rain[(3 * 365 + 1):(5 * 365), 3:54,drop = FALSE]
#' predictions <- predict(fit, data = ensemble_test)
#' idrPit <- pit(predictions, ytest, seed = 123)
#'
#' rawData <- rain[(3 * 365 + 1):(5 * 365), c("HRES", "CTR", paste0("P", 1:50))]
#' rawPit <- pit(rawData, ytest, seed = 123)
#'
#' hist(idrPit, xlab = "Probability Integral Transform",
#'   ylab = "Density", freq = FALSE, main = "Calibrated ensemble")
#' hist(rawPit, xlab = "Probability Integral Transform",
#'   ylab = "Density", freq = FALSE, main = "Raw ensemble")
pit <- function(predictions, y, randomize = TRUE, seed = NULL) {
    UseMethod("pit")
}


#' pit method for class 'idrsd'
#'
#' @method pit idrsd
#' @rdname pit
#' @importFrom stats stepfun
#' @export
pit.idrsd <- function(predictions, y, randomize = TRUE, seed = NULL) {

  # Check input
  if (!is.vector(y, "numeric"))
    stop("'y' must be a numeric vector")
  if (length(y) != length(predictions))
    stop("'y' must have the same length as the predictions")

  pit0 <- function(data, y) {
    # Evaluated CDF (stepfun)
    stats::stepfun(x = data$points, y = c(0, data$cdf))(y)
  }
  pitVals <- mapply(pit0, data = predictions, y = y)

  if (randomize) {
    # Randomization: Find small epsilon to compute left limit of CDF
    sel <- sapply(predictions, nrow) > 1
    if (!any(sel)) {
      eps <- 1
    } else {
      eps <- min(sapply(predictions[sel], function(prd) min(diff(prd$points))))
    }
    lowerPitVals <- mapply(pit0, data = predictions, y = y - eps / 2)
    if (!is.null(seed))
      set.seed(seed)
    sel <- lowerPitVals < pitVals
    if (any(sel)) {
      pitVals[sel] <- stats::runif(sum(sel), min = lowerPitVals[sel],
        max = pitVals[sel])
    }
  }
  pitVals
}

#' pit method for class 'data.frame'
#'
#' @method pit data.frame
#' @rdname pit
#' @importFrom stats ecdf
#' @export
pit.data.frame <- function(predictions, y, randomize = TRUE, seed = NULL) {

  # Check input
  if (!is.vector(y, "numeric"))
    stop("'y' must be a numeric vector")
  n <- nrow(predictions)
  if (length(y) != n)
    stop("'y' must have the same length as the predictions")
  if (!all(sapply(predictions, is.numeric)))
      stop("'predictions' contains non-numeric variables")

  pit0 <- function(data, y) stats::ecdf(data)(y)
  pred <- unname(split(data.matrix(predictions), seq_len(n)))
    # Update split(data.matrix(...)) to asplit
  pitVals <- mapply(pit0, data = pred, y = y)
  if (randomize) {
    # Randomization: Find small epsilon to compute left limit of CDF
    eps <- apply(predictions, 1, stats::dist, method = "manhattan")
    eps <- min(c(eps[eps > 0], 1))
    lowerPitVals <- mapply(pit0, data = pred, y = y - eps / 2)
    if (!is.null(seed))
      set.seed(seed)
    sel <- lowerPitVals < pitVals
    if (any(sel)) {
      pitVals[sel] <- stats::runif(sum(sel), min = lowerPitVals[sel],
        max = pitVals[sel])
    }
  }
  pitVals
}

#' Plot IDRsd predictions
#'
#' @description Plot an IDRsd predictive CDF.
#'
#' @method plot idrsd
#'
#' @param x object of class \code{idrsd} (output of
#'   \code{\link{predict.idrcal}}).
#' @param index index of the prediction in \code{x} for which a plot is desired.
#' @param bounds whether the bounds should be plotted or not (see
#'   \code{\link{predict.idrfit}} for details about the meaning of the bounds).
#' @param col.cdf color of the predictive CDF.
#' @param col.bounds color of the bounds.
#' @param lty.cdf linetype of the predictive CDF.
#' @param lty.bounds linetype of the CDF bounds.
#' @param main main title.
#' @param xlab label for x axis.
#' @param ylab label for y axis.
#' @param ... further arguments to \code{\link{plot.stepfun}} or
#'   \code{\link{plot}}.
#'
#' @return
#' The data based on which the plot is drawn (returned invisible).
#'
#' @seealso \code{\link{predict.idrcal}}
#'
#' @export
#' @importFrom stats stepfun
#' @importFrom graphics plot
#'
#' @examples
#' data("rain")
#' require("graphics")
#'
#'
#' ensemble <- rain[1:(3 * 365), 3:54, drop = FALSE]
#' y <- rain[1:(3 * 365), "obs"]
#'
#' ## Fit IDRsd and plot the predictive CDF
#'
#' fit <- idrsd(y = y, X = ensemble, type = 'ensemble')
#' ensemble_test = rain[(3 * 365 + 1), 3:54,drop = FALSE]
#' pred <- predict(fit, data = ensemble_test)
#' plot(pred)
plot.idrsd <- function(x, index = 1, bounds = TRUE, col.cdf = "black",
   col.bounds = "blue", lty.cdf = 1, lty.bounds = 3, xlab = "Threshold",
   ylab = "CDF", main = "IDR predictive CDF", ...) {
  pred <- x[[index]]
  graphics::plot(stepfun(x = pred$points, y = c(0, pred$cdf)), xlab = xlab,
    ylab = ylab, do.points = FALSE, col = col.cdf, lty = lty.cdf,
    main = main, ...)
  graphics::abline(h = c(0, 1), col = "gray", lty = 3)
  if (bounds && !is.null(pred$upper)) {
    if (any(pred$lower > 0)) {
      graphics::plot(stats::stepfun(x = pred$points, y = c(0, pred$lower)),
        add = TRUE, col = col.bounds, lty = lty.bounds, do.points = FALSE)
    } else {
      graphics::abline(h = 0, col = col.bounds, lty = lty.bounds)
    }

    if (any(pred$upper < 1)) {
      graphics::plot(stats::stepfun(x = pred$points, y = c(0, pred$upper)),
        add = TRUE, col = col.bounds, lty = lty.bounds, do.points = FALSE)
    } else {
      graphics::abline(h = 1, col = col.bounds, lty = lty.bounds)
    }
  }
  invisible(pred)
}

# Uncertainty
crps_unc <- function(obs){
  obs_ensemble <- matrix(rep(obs, length(obs)), nrow = length(obs), ncol = length(obs), byrow = TRUE)
  return(mean(crps_sample(obs, dat = obs_ensemble)))
}

# CRPS for ECDFs
crps_ecdf <- function(y, grid_vals, ecdf) {

  x <- lapply(seq_len(length(y)), function(i) grid_vals)
  p <- lapply(seq_len(nrow(ecdf)), function(i) ecdf[i,])#as.list( t(ecdf) )
  w <- lapply(p, function(x) c(x[1], diff(x)))

  crps0 <- function(y, p, w, x) 2 * sum(w * ((y < x) - p + 0.5 * w) * (x - y))
  mapply(crps0, y = y, p, w = w, x = x)
}

#' CRPS decomposition
#' @description Computes the individual components of the iso-based CRPS Decomposition: MSC, DSC and UNC
#'
#' @param y numeric vector (the response variable).
#' @param X object that fits to specification in \code{type}. For \code{type = 'idr'}, \code{X} must be an IDR object,
#' for \code{type = 'ensemble'}, \code{X} must be a data.frame where columns correspond to ensemble members, for \code{type = 'ecdf'}, \code{X} must be a \code{(n x m)} matrix
#' with \code{n = length(y)} and \code{m = length(grid)}, where each row corresponds to ECDF values evaluated at \code{grid}.
#' For \code{type = 'dis'} \code{X} is empty, for \code{type = 'normal'} and \code{type = 'normal_ab'} \code{X} is a matrix with 2 columns
#' which represent mu and sigma parameters of the normal distribution.
#' @param type default is \code{'ensemble'}. Other possibilities are \code{'idr'}, \code{'ecdf'}, \code{'dis'}, \code{'normal'}, \code{'normal_ab'}
#' @param grid if \code{type == 'ecdf'}, than \code{grid} is vector of threshold values corresponding to ECDF-values in \code{X}
#' @param dis_func if \code{type == 'dis'}, then a cumulative distribution function must be specified along with its distributional parameters.
#' @param pars parameters for quadratic programming optimization (only relevant
#'   if \code{X} has more than one column), set using
#'   \code{\link[osqp]{osqpSettings}}.
#' @param progress display progressbar (\code{TRUE}, \code{FALSE} or \code{1},
#'   \code{0})?
#'
#' @details This function computed the CRPS decomposition of a response vector \emph{y} and a distributional forecast \emph{X}, which can be an ensemble and empirical cumulative distributional function,
#' a normal distribution or any other other closed form CDF specified by its parameters.
#'
#' @return A list of CRPS decomposition components: miscalibration (MSC), discrimination (DSC), uncertainty (UNC) and the original CRPS.
#'
#'
#' @references
#'
#' @export
#' @examples
#'
#'#' data("rain")
#'
#' ## IDRsd based on ensemble forecast
#'
#' ensemble <- rain[1:(3 * 365), 3:54, drop = FALSE]
#' y <- rain[1:(3 * 365), "obs"]
#'
#' ## Compute CRPS decomposition components
#' crps_deco <- isodeco_crps(y = y, X = ensemble, type = 'ensemble')
#' print(crps_deco)
#'
#'
isodeco_crps <- function(y, X=NULL, grid = NULL, dis_func = NULL, type = 'ensemble', # 'dis', 'ecdf', 'normal', 'idr', normal_ab
                      inta = NULL, intb = NULL,
                      pars = osqpSettings(verbose = FALSE, eps_abs = 1e-5,
                                          eps_rel = 1e-5, max_iter = 10000L),
                      progress = TRUE, ...) {

  cali_idr <-  idrsd(y = y, X = X, grid = grid, dis_func = dis_func , type = type, inta = inta, intb = intb, # 'dis', 'ecdf', 'normal', 'idr'
                      pars = pars,progress = progress, org_crps = TRUE, ...)

  cali_preds <- predict(cali_idr)
  cali_crps <- mean(crps(cali_preds, y))

  uncertainty <- crps_unc(y)

  crps_original <- cali_idr$org_crps

  return(list('MCB' = crps_original - cali_crps, 'DSC' =  uncertainty - cali_crps, 'UNC' = uncertainty, 'CRPS' = crps_original))
}



