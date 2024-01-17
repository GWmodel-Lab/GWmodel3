#' Geographically Weighted Summary Statistics
#' 
#' @param formula Regresison model.
#' @param data A `sf` objects.
#' @param bw Bandwidth value
#' @param adaptive Whether the bandwidth value is adaptive or not.
#' @param mode Algorithm mode. Average or correlation.
#' @param quantile Whether to calculate local quantiles.
#' @param kernel Kernel function used.
#' @param longlat Whether the coordinates
#' @param p Power of the Minkowski distance,
#'  default to 2, i.e., Euclidean distance.
#' @param theta Angle in radian to roate the coordinate system, default to 0.
#' @param parallel_method Parallel method.
#' @param parallel_arg Parallel method argument.
#' 
#' @return A `gwssm` object.
#' 
#' @export 
gwss <- function(
    formula,
    data,
    bw,
    adaptive = FALSE,
    mode = c("average", "correlation"),
    quantile = FALSE,
    kernel = c("gaussian", "exp", "bisquare", "tricube", "boxcar"),
    longlat = FALSE,
    p = 2.0,
    theta = 0.0,
    parallel_method = c("no", "omp"),
    parallel_arg = c(0)
) {
    ### Check args
    mode = match.arg(mode)
    kernel = match.arg(kernel)
    parallel_method = match.arg(parallel_method)
    attr(data, "na.action") <- getOption("na.action")

    ### Extract coords
    data <- do.call(na.action(data), args = list(data))
    coords <- as.matrix(sf::st_coordinates(sf::st_centroid(data)))
    if (is.null(coords) || nrow(coords) != nrow(data))
        stop("Missing coordinates.")

    ### Extract variables
    mc <- match.call(expand.dots = FALSE)
    mf <- model.frame(formula = formula(formula), data = sf::st_drop_geometry(data))
    mt <- attr(mf, "terms")
    x <- model.matrix(mt, mf)
    indep_vars <- colnames(x)
    colnames(x) <- indep_vars

    ### Check whether bandwidth is valid.
    if (missing(bw)) {
        stop("Bandwidth is missing.")
    } else if (is.numeric(bw) || is.integer(bw)) {
        if (bw == 0) stop("Bandwidth must be non-zero.")
    } else {
        stop("Bandwidth must be a number.")
    }

    c_results <- gwss_fit(
        x, coords, enum(mode), quantile,
        bw, adaptive, enum(kernel),
        longlat, p, theta,
        enum_list(parallel_method, parallel_types), parallel_arg
    )

    sdf_data <- data.frame()
    if (mode == "average") {
        local_mean <- c_results$mean
        local_var <- c_results$var
        local_sdev <- c_results$sdev
        local_skew <- c_results$skew
        local_cv <- c_results$cv
        colnames(local_mean) <- paste(indep_vars, "_LM")
        colnames(local_var) <- paste(indep_vars, "_LVar")
        colnames(local_sdev) <- paste(indep_vars, "_LSdev")
        colnames(local_skew) <- paste(indep_vars, "_LSkew")
        colnames(local_cv) <- paste(indep_vars, "_LCV")
        sdf_data <- as.data.frame(cbind(local_mean, local_var, local_sdev, local_skew, local_cv))
        if (quantile) {
            local_median = c_results$median
            local_iqr = c_results$iqr
            local_qi = c_results$qi
            colnames(local_median) <- paste(indep_vars, "_LMedian")
            colnames(local_iqr) <- paste(indep_vars, "_IQR")
            colnames(local_qi) <- paste(indep_vars, "_QI")
            sdf_data <- cbind(sdf_data, local_median, local_iqr, local_qi)
        }
    } else {
        local_corr <- c_results$corr
        local_scorr <- c_results$scorr
        var_pairs <- outer(indep_vars, indep_vars, paste, sep = ".")
        var_pairs <- var_pairs[upper.tri(var_pairs)]
        colnames(local_corr) <- paste(var_pairs, "_LCor")
        colnames(local_scorr) <- paste(var_pairs, "_LSCor")
        sdf_data <- cbind(sdf_data, local_median, local_iqr, local_qi)
    }
    sdf_data$geometry <- sf::st_geometry(data)
    sdf <- sf::st_sf(sdf_data)

    ### Return result
    gwssm <- list(
        SDF = sdf,
        args = list(
            x = x,
            coords = coords,
            bw = bw,
            adaptive = adaptive,
            mode = mode,
            quantile = quantile,
            kernel = kernel,
            longlat = longlat,
            p = p,
            theta = theta,
            parallel_method = parallel_method,
            parallel_arg = parallel_arg
        ),
        call = mc,
        indep_vars = indep_vars
    )
    class(gwssm) <- "gwssm"
    gwssm
}

#' Print description of a `gwssm` object
#'
#' @param x An `gwssm` object returned by [gwr_basic()].
#' @param decimal_fmt The format string passing to [base::sprintf()].
#' @inheritDotParams print_table_md
#' 
#' @method print gwssm
#' @rdname print
#' @export
print.gwssm <- function(x, ..., decimal_fmt) {
    if (!inherits(x, "gwssm"))
        stop("It's not a 'gwssm' object.")
    cat("Local summary statistics", fill = T)
    cat("========================", fill = T)
    cat("                 Formula:", deparse(x$call$formula), fill = T)
    cat("                    Data:", deparse(x$call$data), fill = T)
    cat("Number of summary points:", nrow(x$args$x), fill = T)
    cat("         Kernel function:", x$args$kernel, fill = T)
    cat("               Bandwidth:", x$args$bw,
        ifelse(x$args$adaptive, "(Nearest Neighbours)", "(Meters)"),
        ifelse(x$args$optim_bw, paste0(
            "(Optimized accroding to ",
            x$args$optim_bw_criterion,
            ")"
        ), ""), fill = T)
    distance_type <- "Euclidean"
    if (x$args$longlat) distance_type <- "Geodetic"
    else if (x$args$p == 2) distance_type <- "Euclidean"
    else if (x$args$p == 1) distance_type <- "Manhattan"
    else if (is.infinite(x$args$p)) distance_type <- "Chebyshev"
    else distance_type <- "Generalized Minkowski"
    distance_rotated <- (x$args$theta != 0 && x$args$p != 2 && !x$args$longlat)
    cat("                Distance:", distance_type, ifelse(distance_rotated, " (rotated)", ""), fill = T)
    cat("\n", fill = T)
}