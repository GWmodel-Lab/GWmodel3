#' An S4 class to set GWDR configurations
#'
#' @slot bw Bandwidth Value
#'
#' @exportClass GWDRConfig
GWDRConfig <- setClass("GWDRConfig", slots = c(
    bw = "numeric",
    adaptive = "logical",
    kernel = "character"
), prototype = list(
    bw = NA_real_,
    adaptive = FALSE,
    kernel = "gaussian"
))

#' Replicate MGWR config
#'
#' @param x A \linkS4class{GWDRConfig} object.
#' @param \dots Additional arguments.
#' @param times Replication times.
#'
#' @return A list of \linkS4class{GWDRConfig} objects.
#'
#' @examples
#' rep(mgwr_config(36, TRUE, "bisquare"), 4)
#'
#' @name rep-GWDRConfig
NULL

#' @rdname rep-GWDRConfig
#' @export
setMethod(
    "rep",
    signature(x = "GWDRConfig"),
    definition = function(x, ...) {
        mc <- match.call(rep.int)
        mc[[1L]] <- as.name("rep.int")
        eval(mc)
    }
)

#' @rdname rep-GWDRConfig
#' @export
setMethod(
    "rep.int",
    signature(x = "GWDRConfig", times = "numeric"),
    definition = function(x, times = 1) {
        times <- as.integer(floor(times))
        lapply(seq_len(times), function(i) {
            mgwr_config(
                bw = x@bw,
                adaptive = x@adaptive,
                kernel = x@kernel
            )
        })
    }
)

#' Create an instance of GWDRConfig
#'
#' @describeIn GWDRConfig-class
#'
#' @examples
#' mgwr_config(36, TRUE, "bisquare", optim_bw = "AIC")
#'
#' @export
gwdr_config <- function(
    bw = 0.618,
    adaptive = TRUE,
    kernel = c("gaussian", "exp", "bisquare", "tricube", "boxcar")
) {
    kernel <- match.arg(kernel)
    new("GWDRConfig",
        bw = bw,
        adaptive = adaptive,
        kernel = kernel
    )
}
