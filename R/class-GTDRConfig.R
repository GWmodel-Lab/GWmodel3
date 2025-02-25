#' An S4 class to set GTDR configurations
#'
#' @slot bw Bandwidth value.
#' @slot adaptive Whether the bandwidth value is adaptive or not.
#' @slot kernel Kernel function used.
#' 
#' @exportClass GTDRConfig
GTDRConfig <- setClass("GTDRConfig", slots = c(
    bw = "numeric",
    adaptive = "logical",
    kernel = "character"
), prototype = list(
    bw = NA_real_,
    adaptive = FALSE,
    kernel = "gaussian"
))

#' Replicate GTDR config
#'
#' @param x A \linkS4class{GTDRConfig} object.
#' @param \dots Additional arguments.
#' @param times Replication times.
#'
#' @return A list of \linkS4class{GTDRConfig} objects.
#'
#' @examples
#' rep(gtdr_config(36, TRUE, "bisquare"), 4)
#'
#' @name rep-GTDRConfig
NULL

#' @rdname rep-GTDRConfig
#' @export
setMethod(
    "rep",
    signature(x = "GTDRConfig"),
    definition = function(x, ...) {
        mc <- match.call(rep.int)
        mc[[1L]] <- as.name("rep.int")
        eval(mc)
    }
)

#' @rdname rep-GTDRConfig
#' @export
setMethod(
    "rep.int",
    signature(x = "GTDRConfig", times = "numeric"),
    definition = function(x, times = 1) {
        times <- as.integer(floor(times))
        lapply(seq_len(times), function(i) {
            gtdr_config(
                bw = x@bw,
                adaptive = x@adaptive,
                kernel = x@kernel
            )
        })
    }
)

#' Create an instance of GTDRConfig.
#'
#' @param bw Bandwidth value.
#' @param adaptive Whether the bandwidth value is adaptive or not.
#' @param kernel Kernel function used.
#'
#' @examples
#' gtdr_config(36, TRUE, "bisquare")
#' 
#' @importFrom methods new
#'
#' @export
gtdr_config <- function(
    bw = 0.618,
    adaptive = TRUE,
    kernel = c("gaussian", "exp", "bisquare", "tricube", "boxcar")
) {
    kernel <- match.arg(kernel)
    new("GTDRConfig",
        bw = bw,
        adaptive = adaptive,
        kernel = kernel
    )
}

gtdr_bw_criterion_enums <- c(
    "CV",
    "AIC"
)
