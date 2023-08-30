#' CRPS decompositiom
#'
#' @section How does it work?:
#' @useDynLib isodisregAFSD, .registration = TRUE
NULL


#' @keywords internal
.onUnload <- function (libpath) {
  library.dynam.unload("isodisregAFSD", libpath)
}
