#' Calibrate a basic GWR model
#' 
#' @param formula Regresison model.
#' @param data A `sf` objects.
#' @param bw Bandwidth value. Depends on the value of `adaptive`. 
#' @param adaptive Whether the bandwidth value is adaptive or not.
#' @param kernel Kernel function used.
#' @param longlat Whether the coordinates
#' @param p Power of the Minkowski distance, default is 2, i.e. the Euclidean distance.
#' @param theta Angle in radian to roate the coordinate system, default to 0.
#' @param hatmatrix If TRUE, great circle will be caculated
#' @param parallel_method Parallel method.
#' @param parallel_arg Parallel method argument.
#' 
#' @return A `gwrm` object.
#' 
#' @examples 
#' data(LondonHP)
#' gwr_basic(PURCHASE~FLOORSZ+UNEMPLOY, LondonHP, 64, TRUE)
#' 
#' @export
gwr_basic <- function(
    formula,
    data,
    bw,
    adaptive = FALSE,
    kernel = c("gaussian", "exp", "bisquare", "tricube", "boxcar"),
    longlat = FALSE,
    p = 2.0,
    theta = 0.0,
    hatmatrix = TRUE,
    parallel_method = c("none", "omp", "cuda", "cluster"),
    parallel_arg = c(0)) 
{
    ### Check args
    kernel = match.arg(kernel)
    parallel_method = match.arg(parallel_method)

    ### Extract coords
    coords <- as.matrix(sf::st_coordinates(sf::st_centroid(data)))
    if (is.null(coords) || nrow(coords) != nrow(data))
        stop("Missing coordinates.")
    
    ### Extract variables
    fields <- all.vars(formula)
    xfields <- fields[-1]
    yfields <- fields[1]
    x <- as.matrix(sf::st_drop_geometry(data[, xfields]))
    y <- as.matrix(sf::st_drop_geometry(data[, yfields]))

    ### Call solver
    c_result <- .c_gwr_basic(
        x, y, coords, xfields, yfields,
        bw, adaptive, kernel, longlat, p, theta,
        hatmatrix, parallel_method, parallel_arg
    )

    ### Return result
    class(c_result) <- "gwrm"
    c_result
}
