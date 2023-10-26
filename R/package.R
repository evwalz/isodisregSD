#' Isotonic distributional Regression (IDR) for distributional data
#'
#' Isotonic distributional Regression (IDR) for distributional data is a nonparametric method to
#' estimate conditional distributions under stochastic order constraints.
#'
#' @section How does it work?:
#' Link to Preprint on ArXiv as soon as available
#' @docType package
#' @name isodisregSD-package
#' @useDynLib isodisregSD, .registration = TRUE
NULL


#' @keywords internal
.onUnload <- function (libpath) {
  library.dynam.unload("isodisregSD", libpath)
}
