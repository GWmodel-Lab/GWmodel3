#' Geographically Weighted Average
#' 
#' @param data A `sf` objects.
#' @param vars a vector of variable names to be summarized.
#' @param bw Bandwidth value.
#' @param adaptive Whether the bandwidth value is adaptive or not.
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
#' @examples
#' data(LondonHP)
#' m <- gw.average(LondonHP, c("PURCHASE","FLOORSZ","UNEMPLOY"), 64, TRUE)
#' gw.average(LondonHP, "PURCHASE,FLOORSZ,UNEMPLOY", 64, TRUE)
#' 
#' @export 
gw.average <- function(
    data,
    vars,
    bw,
    adaptive = FALSE,
    quantile = FALSE,
    kernel = c("gaussian", "exp", "bisquare", "tricube", "boxcar"),
    longlat = FALSE,
    p = 2.0,
    theta = 0.0,
    parallel_method = c("no", "omp"),
    parallel_arg = c(0)
) {
    ### Check args
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

    if (inherits(data, "sf")) {
        df <- sf::st_drop_geometry(data)
    } else {
        stop("Input data should be an sf object.")
    }

    if (is.character(vars)) {
        if (length(vars) == 1 && grepl(",", vars)) {
            vars <- unlist(strsplit(vars, ","))
        }
        # only one variable
    } else if (is.vector(vars) && !is.character(vars)) {
        stop("Input 'vars' should be either a character vector or a comma-separated string.")
    }

    len.var <- length(vars)
    col.nm <- colnames(df)
    var.idx <- match(vars, col.nm)[!is.na(match(vars, col.nm))]
    if (length(var.idx) == 0) 
        stop("All variables input doesn't match with data")
    x <- df[, var.idx]
    x <- as.matrix(x)
    indep_vars <- colnames(x)
    # colnames(x) <- indep_vars

    ### Check whether bandwidth is valid.
    if (missing(bw)) {
        stop("Bandwidth is missing.")
    } else if (is.numeric(bw) || is.integer(bw)) {
        if (bw == 0) stop("Bandwidth must be non-zero.")
    } else {
        stop("Bandwidth must be a number.")
    }

    c_results <- gw_average(
        x, coords, quantile,
        bw, adaptive, enum(kernel),
        longlat, p, theta,
        enum_list(parallel_method, parallel_types), parallel_arg
    )

    sdf_data <- data.frame()
    ## copy results
    {
        local_mean <- c_results$mean
        local_var <- c_results$var
        local_sdev <- c_results$sdev
        local_skew <- c_results$skew
        local_cv <- c_results$cv
        colnames(local_mean) <- paste0(indep_vars, "_LM")
        colnames(local_var) <- paste0(indep_vars, "_LVar")
        colnames(local_sdev) <- paste0(indep_vars, "_LSdev")
        colnames(local_skew) <- paste0(indep_vars, "_LSkew")
        colnames(local_cv) <- paste0(indep_vars, "_LCV")
        sdf_data <- as.data.frame(cbind(local_mean, local_var, local_sdev, local_skew, local_cv))
        if (quantile) {
            local_median = c_results$median
            local_iqr = c_results$iqr
            local_qi = c_results$qi
            colnames(local_median) <- paste0(indep_vars, "_LMedian")
            colnames(local_iqr) <- paste0(indep_vars, "_IQR")
            colnames(local_qi) <- paste0(indep_vars, "_QI")
            sdf_data <- cbind(sdf_data, local_median, local_iqr, local_qi)
        }
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
            mode = "Average",
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
#' @importFrom sf st_drop_geometry
#' @method print gwssm
#' @rdname print
#' @export
print.gwssm <- function(x, ..., decimal_fmt) {
    if (!inherits(x, "gwssm"))
        stop("It's not a 'gwssm' object.")
    cat("Local summary statistics", fill = T)
    cat("========================", fill = T)
    # cat("                 Formula:", deparse(x$call$formula), fill = T)
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

    cat("SDF Summary", fill = T)
    cat("===========", fill = T)
    res <- st_drop_geometry(x$SDF)
    print(as.data.frame(t(sapply(res, summary))))
}

#' @describeIn gwss Plot the result of basic GWSS model.
#' 
#' @param x A "gwssm" object
#' @param y Ignored
#' @param columnsColumn names to plot.
#'  If it is missing or non-character value, all coefficient columns are plotted.
#' @param \dots Additional arguments passing to [sf::plot()].
#' @method plot gwssm
#' 
#' @examples 
#' plot(m)
#' 
#' @export
plot.gwssm <- function(x, y, ..., columns) {
    if (!inherits(x, "gwssm"))
        stop("It's not a 'gwssm' object.")
    
    sdf <- sf::st_as_sf(x$SDF)
    sdf_colnames <- names(sf::st_drop_geometry(x$SDF))
    if (!missing(columns) && is.character(columns)) {
        valid_columns <- intersect(columns, sdf_colnames)
        if (length(valid_columns) > 0) {
            sdf <- sdf[valid_columns]
        }
    }
    plot(sdf, ...)
}
