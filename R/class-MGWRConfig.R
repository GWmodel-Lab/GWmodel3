#' An S4 class to set Multiscale GWR configurations
#'
#' @slot bw Bandwidth Value
#'
#' @exportClass MGWRConfig
MGWRConfig <- setClass("MGWRConfig", slots = c(
    bw = "numeric",
    adaptive = "logical",
    kernel = "character",
    longlat = "logical",
    p = "numeric",
    theta = "numeric",
    centered = "logical",
    optim_bw = "character",
    optim_threshold = "numeric"
), prototype = list(
    bw = NA_real_,
    adaptive = FALSE,
    kernel = "gaussian",
    longlat = FALSE,
    p = 2.0,
    theta = 0.0,
    centered = FALSE,
    optim_bw = "AIC",
    optim_threshold = 1e-5
))

setValidity("MGWRConfig", function(object) {
    if (is.na(object@bw) && object@optim_bw == "no") {
        "@bw cannot be NA when @optim_bw is 'no'"
    } else {
        TRUE
    }
})

#' Replicate MGWR config
#'
#' @param x A \linkS4class{MGWRConfig} object.
#' @param \dots Additional arguments.
#' @param times Replication times.
#'
#' @return A list of \linkS4class{MGWRConfig} objects.
#'
#' @examples
#' rep(mgwr_config(36, TRUE, "bisquare"), 4)
#'
#' @name rep-MGWRConfig
NULL

#' @rdname rep-MGWRConfig
#' @export
setMethod(
    "rep",
    signature(x = "MGWRConfig"),
    definition = function(x, ...) {
        mc <- match.call(rep.int)
        mc[[1L]] <- as.name("rep.int")
        eval(mc)
    }
)

#' @rdname rep-MGWRConfig
#' @export
setMethod(
    "rep.int",
    signature(x = "MGWRConfig", times = "numeric"),
    definition = function(x, times = 1) {
        times <- as.integer(floor(times))
        lapply(seq_len(times), function(i) {
            mgwr_config(
                bw = x@bw,
                adaptive = x@adaptive,
                kernel = x@kernel,
                longlat = x@longlat,
                p = x@p,
                theta = x@theta,
                centered = x@centered,
                optim_bw = x@optim_bw,
                optim_threshold = x@optim_threshold
            )
        })
    }
)

#' Create an instance of MGWRConfig
#'
#' @describeIn MGWRConfig-class
#'
#' @examples
#' mgwr_config(36, TRUE, "bisquare", optim_bw = "AIC")
#'
#' @export
mgwr_config <- function(
    bw = NA_real_,
    adaptive = FALSE,
    kernel = c("gaussian", "exp", "bisquare", "tricube", "boxcar"),
    longlat = FALSE,
    p = 2.0,
    theta = 0.0,
    centered = FALSE,
    optim_bw = c("AIC", "CV", "no"),
    optim_threshold = 1e-5
) {
    kernel <- match.arg(kernel)
    optim_bw <- match.arg(optim_bw)
    if (is.na(bw) && optim_bw == "no") {
        stop("Cannot specify a NA value as specified bandwidth!")
    }

    new("MGWRConfig",
        bw = bw,
        adaptive = adaptive,
        kernel = kernel,
        longlat = longlat,
        p = p,
        theta = theta,
        centered = centered,
        optim_bw = optim_bw,
        optim_threshold = optim_threshold
    )
}
