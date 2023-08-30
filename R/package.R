#' Isotonic distributional Regression (IDR) for distributional data
#'
#' Isotonic distributional Regression (IDR) is a nonparametric method to
#' estimate conditional distributions under monotonicity constraints.
#'
#' @section How does it work?:
#' Provide some more detailes
#' @docType package
#' @name isodisregAFSD-package
#' @useDynLib isodisregAFSD, .registration = TRUE
NULL


#' @keywords internal
.onUnload <- function (libpath) {
  library.dynam.unload("isodisregAFSD", libpath)
}
