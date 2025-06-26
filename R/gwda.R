#' Geographically Weighted Discriminant Analysis
#'
#' @param formula Regresison model.
#' @param data A `sf` objects.
#' @param bw Bandwidth value
#' @param adaptive Whether the bandwidth value is adaptive or not.
#' @param quantile Whether to calculate local quantiles.
#' @param kernel Kernel function used.
#' @param longlat Whether the coordinates
#' @param p Power of the Minkowski distance,
#'  default to 2, i.e., Euclidean distance.
#' @param theta Angle in radian to roate the coordinate system, default to 0.
#' @param method Discriminant method for GWDA, can be set as
#'  - `wqda`, weighted quadratic discriminant analysis
#'  - `wlda`, weighted linear discriminant analysis
#'  If not provided, this parameter will be set to `wqda`.
#' @param parallel_method Parallel method.
#' @param parallel_arg Parallel method argument.
#'
#' @return A `gwdam` object.
#'
#' @export
gwda <- function(
    formula,
    data,
    bw,
    adaptive = FALSE,
    kernel = c("gaussian", "exp", "bisquare", "tricube", "boxcar"),
    longlat = FALSE,
    p = 2.0,
    theta = 0.0,
    method = c("wqda", "wlda"),
    parallel_method = c("no", "omp"),
    parallel_arg = c(0))
{
    ### Check args
    kernel <- match.arg(kernel)
    parallel_method <- match.arg(parallel_method)
    attr(data, "na.action") <- getOption("na.action")

    ### Extract coords
    data <- do.call(na.action(data), args = list(data))
    coords <- as.matrix(sf::st_coordinates(sf::st_centroid(data)))
    if (is.null(coords) || nrow(coords) != nrow(data)) {
        stop("Missing coordinates.")
    }

    ### Extract variables
    mc <- match.call(expand.dots = FALSE)
    mf <- model.frame(formula = formula(update(formula, ~ . - 1)), data = sf::st_drop_geometry(data))
    mt <- attr(mf, "terms")
    y <- model.extract(mf, "response")
    y <- as.character(y)
    x <- model.matrix(mt, mf)
    dep_var <- as.character(attr(terms(mf), "variables")[[2]])
    indep_vars <- colnames(x)

    ### Check whether bandwidth is valid.
    if (missing(bw)) {
        stop("Bandwidth is missing.")
    } else if (is.numeric(bw) || is.integer(bw)) {
        if (bw == 0) stop("Bandwidth must be non-zero.")
    } else {
        stop("Bandwidth must be a number.")
    }

    method <- match.arg(method)
    wqda_method <- ifelse(method == "wqda", TRUE, FALSE)

    c_results <- gwda_cal(
        x, y, coords, bw,
        adaptive, enum(kernel),
        longlat, p, theta,
        wqda_method,
        enum_list(parallel_method, parallel_types), parallel_arg
    )

    sdf_data <- data.frame()
    res <- c_results$res
    group <- as.data.frame(c_results$group)
    probs <- c_results$probs
    pmax <- c_results$pmax
    entropy <- c_results$entropy

    colnames(res) <- indep_vars
    colnames(group) <- "groups"
    colnames(probs) <- paste0(indep_vars, "_probs")
    colnames(pmax) <- "pmax"
    colnames(entropy) <- "entropy"
    sdf_data <- as.data.frame(cbind(res, group, probs, pmax, entropy))
    
    sdf_data$geometry <- sf::st_geometry(data)
    sdf <- sf::st_sf(sdf_data)

    ### Return result
    gwdam <- list(
        SDF = sdf,
        args = list(
            x = x,
            y = y,
            correctRate = c_results$correctRate,
            coords = coords,
            bw = bw,
            adaptive = adaptive,
            kernel = kernel,
            longlat = longlat,
            p = p,
            theta = theta,
            method = method,
            parallel_method = parallel_method,
            parallel_arg = parallel_arg
        ),
        call = mc,
        indep_vars = indep_vars,
        dep_var = dep_var
    )
    class(gwdam) <- "gwdam"
    gwdam
}

#' Print description of a `gwdam` object
#'
#' @param x An `gwdam` object returned by [gwr_basic()].
#' @param decimal_fmt The format string passing to [base::sprintf()].
#' @inheritDotParams print_table_md
#'
#' @method print gwdam
#' @rdname print
#' @export
print.gwdam <- function(x, ..., decimal_fmt) {
    if (!inherits(x, "gwdam")) {
        stop("It's not a 'gwdam' object.")
    }

    cat("   ***********************************************************************\n")
    cat("   *                         Package   GWmodel3                          *\n")
    cat("   ***********************************************************************\n")
    cat("   *            Results of Geographically Discriminant Analysis          *\n")
    cat("   ***********************************************************************\n")
    cat("\n   *********************Model Calibration Information*********************\n")

    cat("   Formula:", deparse(x$call$formula), fill = T)
    cat("   Data:", deparse(x$call$data), fill = T)
    cat("   Method:", x$args$method, fill = T)
    cat("   Number of summary points:", nrow(x$args$x), fill = T)
    cat("   Kernel function:", x$args$kernel, fill = T)
    cat("   Bandwidth:", x$args$bw,
        ifelse(x$args$adaptive, "(Nearest Neighbours)", "(Meters)"),
        fill = T
    )
    distance_type <- "Euclidean"
    if (x$args$longlat) {
        distance_type <- "Geodetic"
    } else if (x$args$p == 2) {
        distance_type <- "Euclidean"
    } else if (x$args$p == 1) {
        distance_type <- "Manhattan"
    } else if (is.infinite(x$args$p)) {
        distance_type <- "Chebyshev"
    } else {
        distance_type <- "Generalized Minkowski"
    }
    distance_rotated <- (x$args$theta != 0 && x$args$p != 2 && !x$args$longlat)
    cat("   Distance:", distance_type, ifelse(distance_rotated, " (rotated)", ""), fill = T)
    res <- st_drop_geometry(x$SDF)
    group_counts <- as.data.frame(table(res$groups))
    cat("   Discrimination Result:", fill = T)
    colnames(group_counts) <- c("Group", "Count")
    output <- capture.output(print(group_counts, row.names = FALSE))
    output <- paste0("  ", output)
    cat(output, sep = "\n")

    cat("\n   ***********************Local Summary Statistics************************", fill = T)
    numeric_res <- res[sapply(res, is.numeric)]
    df <- as.data.frame(t(sapply(numeric_res, summary)))
    output <- capture.output(print(df))
    output <- paste0("   ", output)
    cat(output, sep = "\n")
    cat("   ***********************************************************************\n")
}