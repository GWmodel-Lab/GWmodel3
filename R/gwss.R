#' Geographically Weighted Average
#' 
#' @param data A `sf` objects.
#' @param vars a vector of variable names to be summarized.
#' @param bw Bandwidth value, need to be specified.
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
    colnames(x) <- indep_vars

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
            vars1 = x,
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
        name_vars1 = indep_vars
    )
    class(gwssm) <- "gwssm"
    gwssm
}

#' Geographically Weighted Correlation
#' 
#' @param data A `sf` objects.
#' @param vars1 a vector of variable names to be calculated.
#' @param vars2 a vector of variable names as the responsed.
#' @param bws Bandwidth value list, auto-select if not provided.
#' @param adaptive Whether the bandwidth value is adaptive or not.
#' @param kernel Kernel function used.
#' @param approach Bandwidth selection.
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
#' m <- gw.correlation(LondonHP, c("PURCHASE","PROF"),c("FLOORSZ","UNEMPLOY"),kernel="bisquare", adaptive=TRUE)
#' gw.correlation(LondonHP, c("PURCHASE"),c("FLOORSZ","UNEMPLOY"),c(50,50), adaptive=TRUE)
#' 
#' @export 
gw.correlation <- function(
    data,
    vars1,
    vars2,
    bws = NULL,
    adaptive = FALSE,
    kernel = c("gaussian", "exp", "bisquare", "tricube", "boxcar"),
    approach = c("CV", "AIC"),
    longlat = FALSE,
    p = 2.0,
    theta = 0.0,
    parallel_method = c("no", "omp"),
    parallel_arg = c(0)
) {
    ### Check args
    kernels = match.arg(kernel)
    parallel_method = match.arg(parallel_method)
    approach = match.arg(approach)
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

    if (!is.character(vars1) || length(vars1) == 0 ) {
        stop("Input 'vars1' should be a character vector")
    }
    if (!is.character(vars2) || length(vars2) == 0 ) {
        stop("Input 'vars2' should be a character vector")
    }

    len.var <- length(vars1) * length(vars2)
    col.nm <- colnames(df)

    var.idx1 <- match(vars1, col.nm)[!is.na(match(vars1, col.nm))]
    if (length(var.idx1) == 0) 
        stop("All variables input doesn't match with data")
    x <- df[, var.idx1]
    x <- as.matrix(x)
    colnames(x) <- vars1
    vars1 <- colnames(x)
    var.idx2 <- match(vars2, col.nm)[!is.na(match(vars2, col.nm))]
    if (length(var.idx2) == 0) 
        stop("All variables input doesn't match with data")
    y <- df[, var.idx2]
    y <- as.matrix(y)
    colnames(y) <- vars2
    vars2 <- colnames(y)

    ### Check bandwidth.
    initial_type <- c(rep("Specified", length(bws)), rep("Null", ifelse(len.var>=length(bws),len.var - length(bws),len.var)))
    bw_select_type <- c(rep("Specified", length(bws)), rep("AutoSelect", ifelse(len.var>=length(bws),len.var - length(bws),len.var)))
    if (missing(approach))  approach <- "AIC"
    optim_bw_criterion <- rep(match.arg(approach), len.var)

    if (missing(bws)) {
        bws <- c(rep(1, len.var))
        # optim_bw <- TRUE
    } else {
        if (is.numeric(bws) || is.integer(bws)) {
            # optim_bw <- FALSE
            if (any(bws <= 0)) stop("Bandwidth must be positive numbers.")
        } else {
            stop("Bandwidth must be a number.")
        }
    }

    if (length(bws) > len.var) {
        bws <- bws[1:len.var]
        # optim_bw <- FALSE
    }

    ### Config kernels
    vkernel <- c(rep(kernel, len.var))
    # vkernel <- c(rep(which(kernel == kernels) - 1,len.var))
    vadaptive <- c(rep(adaptive, len.var))
    vlonglat <- c(rep(longlat, len.var))
    vp <- c(rep(p, len.var))
    vtheta <- c(rep(theta, len.var))

    c_results <- tryCatch(gw_correlation(
        x, y, coords,
        bws, vadaptive, enum(vkernel, kernel_enums),
        vlonglat, vp, vtheta,
        enum(initial_type,initial_enums), enum(optim_bw_criterion,criterion_enums),
        enum_list(parallel_method, parallel_types), parallel_arg
    ), error = function (e) {
        stop("Error:", conditionMessage(e))
    })

    sdf_data <- data.frame()
    ## copy results
    {
        local_cov <- c_results$Cov
        local_corr <- c_results$Corr
        local_scorr <- c_results$SCorr
        colnames(local_cov) <- paste(paste(rep(vars1, each = length(vars2)), rep(vars2, times = length(vars1)), sep = "_"), "LCov", sep = "_")
        colnames(local_corr) <- paste(paste(rep(vars1, each = length(vars2)), rep(vars2, times = length(vars1)), sep = "_"), "LCorr", sep = "_")
        colnames(local_scorr) <- paste(paste(rep(vars1, each = length(vars2)), rep(vars2, times = length(vars1)), sep = "_"), "LSCorr", sep = "_")
        sdf_data <- as.data.frame(cbind(local_cov, local_corr, local_scorr))
    } 
    sdf_data$geometry <- sf::st_geometry(data)
    sdf <- sf::st_sf(sdf_data)

    ### Return result
    gwssm <- list(
        SDF = sdf,
        args = list(
            vars1 = x,
            vars2 = y,
            coords = coords,
            bw = c_results$bw_value,
            mode = "Correlation",
            adaptive = adaptive,
            kernel = kernels,
            longlat = longlat,
            p = p,
            theta = theta,
            # optim_bw = optim_bw,
            bw_configuration = bw_select_type,
            optim_bw_criterion = approach,
            parallel_method = parallel_method,
            parallel_arg = parallel_arg
        ),
        call = mc,
        name_vars1 = vars1,
        name_vars2 = vars2
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
    cat("   ***********************************************************************\n")
    cat("   *                       Package   GWmodel                             *\n")
    cat("   ***********************************************************************\n")
    cat("   ***********************Calibration information*************************\n")
    cat("                   Local summary statistics of GW", fill = T, x$args$mode)
    # cat("   =======================================================================\n", fill = T)
    cat("   Data:", deparse(x$call$data), fill = T)
    cat("   Number of summary points:", nrow(x$args$vars1), fill = T)
    cat("   Kernel function:", x$args$kernel, fill = T)
    switch(x$args$mode,
            "Average" = {
                cat("   Bandwidth:", x$args$bw,
                    ifelse(x$args$adaptive, "(Nearest Neighbours)", "(Meters)"),"\n")
            },
            "Correlation" = {
                var_width <- max(nchar(x$name_vars2))
                bw_matrix <- matrix(x$args$bw, nrow = length(x$name_vars1), ncol = length(x$name_vars2), byrow = TRUE)
                cat("   Bandwidth:\n")
                for (i in 1:length(x$name_vars1)) {
                    cat("             ", 
                        sprintf(x$name_vars1[i]),
                        paste(sprintf("%-*s", var_width, x$name_vars2), collapse = " "),fill = TRUE)
                    cat("             ", 
                        paste(sprintf("%-*s", var_width, "")),
                        paste(sprintf("%-*s", var_width, bw_matrix[i, ]), collapse = ""),  
                        ifelse(x$args$adaptive, "(Nearest Neighbours)", "(Meters)"),
                        # ifelse(x$args$optim_bw, paste0("(Optimized according to ", x$args$optim_bw_criterion, ")"), ""),
                        fill = TRUE
                    )
                }
            }
    )
    distance_type <- "Euclidean"
    if (x$args$longlat) distance_type <- "Geodetic"
    else if (x$args$p == 2) distance_type <- "Euclidean"
    else if (x$args$p == 1) distance_type <- "Manhattan"
    else if (is.infinite(x$args$p)) distance_type <- "Chebyshev"
    else distance_type <- "Generalized Minkowski"
    distance_rotated <- (x$args$theta != 0 && x$args$p != 2 && !x$args$longlat)
    cat("   Distance metric:", distance_type, ifelse(distance_rotated, " (rotated)", "distance metric is used."), fill = T)
    
    cat("\n   ************************Local Summary Statistics:**********************\n", fill = T)
    # cat("   =======================================================================\n", fill = T)
    res <- st_drop_geometry(x$SDF)
    print(as.data.frame(t(sapply(res, summary))))
    cat("   ***********************************************************************\n")
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


criterion_enums <- c(
    "CV",
    "AIC"
)

initial_enums <- c(
    "Null",
    "Specified"
)