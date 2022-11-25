#' Generic method for model auto-selection.
#' @param x Geographically Weighted Model
#' @param \dots Other arguments passing to implementation.
#' 
#' @export
model_sel <- function(x, ...) {
    UseMethod("model_sel")
}
