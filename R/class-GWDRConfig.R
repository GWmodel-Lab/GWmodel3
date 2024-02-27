#' An S4 class to set SDR configurations
#'
#' @slot bw Bandwidth value.
#' @slot adaptive Whether the bandwidth value is adaptive or not.
#' @slot kernel Kernel function used.
#' 
#' @exportClass SDRConfig
SDRConfig <- setClass("SDRConfig", slots = c(
    bw = "numeric",
    adaptive = "logical",
    kernel = "character"
), prototype = list(
    bw = NA_real_,
    adaptive = FALSE,
    kernel = "gaussian"
))

#' Replicate SDR config
#'
#' @param x A \linkS4class{SDRConfig} object.
#' @param \dots Additional arguments.
#' @param times Replication times.
#'
#' @return A list of \linkS4class{SDRConfig} objects.
#'
#' @examples
#' rep(sdr_config(36, TRUE, "bisquare"), 4)
#'
#' @name rep-SDRConfig
NULL

#' @rdname rep-SDRConfig
#' @export
setMethod(
    "rep",
    signature(x = "SDRConfig"),
    definition = function(x, ...) {
        mc <- match.call(rep.int)
        mc[[1L]] <- as.name("rep.int")
        eval(mc)
    }
)

#' @rdname rep-SDRConfig
#' @export
setMethod(
    "rep.int",
    signature(x = "SDRConfig", times = "numeric"),
    definition = function(x, times = 1) {
        times <- as.integer(floor(times))
        lapply(seq_len(times), function(i) {
            sdr_config(
                bw = x@bw,
                adaptive = x@adaptive,
                kernel = x@kernel
            )
        })
    }
)

#' Create an instance of SDRConfig.
#'
#' @param bw Bandwidth value.
#' @param adaptive Whether the bandwidth value is adaptive or not.
#' @param kernel Kernel function used.
#'
#' @examples
#' sdr_config(36, TRUE, "bisquare")
#' 
#' @importFrom methods new
#'
#' @export
sdr_config <- function(
    bw = 0.618,
    adaptive = TRUE,
    kernel = c("gaussian", "exp", "bisquare", "tricube", "boxcar")
) {
    kernel <- match.arg(kernel)
    new("SDRConfig",
        bw = bw,
        adaptive = adaptive,
        kernel = kernel
    )
}

sdr_bw_criterion_enums <- c(
    "CV",
    "AIC"
)
