#' An S4 class to set GW Correlation configurations
#'
#' @slot bw Bandwidth value.
#' @slot adaptive Whether the bandwidth value is adaptive or not.
#' @slot kernel Kernel function used.
#' @slot longlat Whether the coordinates.
#' @slot p Power of the Minkowski distance,
#'  default to 2, i.e., Euclidean distance.
#' @slot theta Angle in radian to roate the coordinate system, default to 0.
#' @slot optim_bw Whether optimize bandwidth after selecting models.
#'  Avaliable values are `no`, `AIC`, and `CV`.
#'  If `no` is specified, the bandwidth specified by argument `bw`
#'  is used in calibrating selected models.
#'
#' @exportClass GWCorrConfig
GWCorrConfig <- setClass("GWCorrConfig", slots = c(
    bw = "numeric",
    adaptive = "logical",
    kernel = "character",
    longlat = "logical",
    p = "numeric",
    theta = "numeric",
    optim_bw = "character"
), prototype = list(
    bw = NA_real_,
    adaptive = FALSE,
    kernel = "gaussian",
    longlat = FALSE,
    p = 2.0,
    theta = 0.0,
    optim_bw = "CV"
))

setValidity("GWCorrConfig", function(object) {
    if (is.na(object@bw) && object@optim_bw == "no") {
        "@bw cannot be NA when @optim_bw is 'no'"
    } else {
        TRUE
    }
})

#' Replicate MGWR config
#'
#' @param x A \linkS4class{GWCorrConfig} object.
#' @param \dots Additional arguments.
#' @param times Replication times.
#'
#' @return A list of \linkS4class{GWCorrConfig} objects.
#'
#' @examples
#' rep(gwcorr_config(36, TRUE, "bisquare"), 4)
#'
#' @name rep-GWCorrConfig
NULL

#' @rdname rep-GWCorrConfig
#' @export
setMethod(
    "rep",
    signature(x = "GWCorrConfig"),
    definition = function(x, ...) {
        mc <- match.call(rep.int)
        mc[[1L]] <- as.name("rep.int")
        eval(mc)
    }
)

#' @rdname rep-GWCorrConfig
#' @export
setMethod(
    "rep.int",
    signature(x = "GWCorrConfig", times = "numeric"),
    definition = function(x, times = 1) {
        times <- as.integer(floor(times))
        lapply(seq_len(times), function(i) {
            gwcorr_config(
                bw = x@bw,
                adaptive = x@adaptive,
                kernel = x@kernel,
                longlat = x@longlat,
                p = x@p,
                theta = x@theta,
                optim_bw = x@optim_bw
            )
        })
    }
)

#' Create an instance of GWCorrConfig
#' 
#' @param bw Bandwidth value.
#' @param adaptive Whether the bandwidth value is adaptive or not.
#' @param kernel Kernel function used.
#' @param longlat Whether the coordinates.
#' @param p Power of the Minkowski distance,
#'  default to 2, i.e., Euclidean distance.
#' @param theta Angle in radian to roate the coordinate system, default to 0.
#' @param optim_bw Whether optimize bandwidth after selecting models.
#'  Avaliable values are `no`, `AIC`, and `CV`.
#'  If `no` is specified, the bandwidth specified by argument `bw`
#'  is used in calibrating selected models.
#'
#' @examples
#' gwcorr_config(36, TRUE, "bisquare", optim_bw = "AIC")
#' 
#' @importFrom methods new
#'
#' @export
#' 
gwcorr_config <- function(
    bw = NA_real_,
    adaptive = FALSE,
    kernel = c("gaussian", "exp", "bisquare", "tricube", "boxcar"),
    longlat = FALSE,
    p = 2.0,
    theta = 0.0,
    optim_bw = c("CV", "AIC", "no")
) {
    kernel <- match.arg(kernel)
    optim_bw <- match.arg(optim_bw)
    if (is.na(bw) && optim_bw == "no") {
        stop("Cannot specify a NA value as specified bandwidth!")
    }

    new("GWCorrConfig",
        bw = bw,
        adaptive = adaptive,
        kernel = kernel,
        longlat = longlat,
        p = p,
        theta = theta,
        optim_bw = optim_bw
    )
}

gwcorr_bw_criterion_enums <- c(
    "CV",
    "AIC"
)

gwcorr_initial_enums <- c(
    "Null",
    "Initial",
    "Specified"
)
