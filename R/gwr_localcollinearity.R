#' Calibrate a basic GWR model
#'
#' @param formula Regresison model.
#' @param data A `sf` objects.
#' @param bw Either a value to set the size of bandwidth,
#'  or one of the following characters to set the criterion for
#'  bandwidth auto-optimization process.
#'  - `AIC`
#'  - `CV`
#'  Note that if `NA` or other non-numeric value is setted,
#'  this parameter will be reset to `Inf`.
#' @param adaptive Whether the bandwidth value is adaptive or not.
#' @param kernel Kernel function used.
#' @param longlat Whether the coordinates
#' @param p Power of the Minkowski distance,
#'  default to 2, i.e., Euclidean distance.
#' @param theta Angle in radian to roate the coordinate system, default to 0.
#' @param optim_bw_range Bounds on bandwidth optimization, a vector of two numeric elements.
#'  Set to `NA_real_` to enable default values selected by the algorithm.
#' @param hatmatrix If TRUE, great circle will be caculated.
#' @param parallel_method Parallel method.
#' @param parallel_arg Parallel method argument.
#' @param verbose Whether to print additional information.
#'
#' @return A `gwlcrm` object.
#' 
#' @details
#' ## Parallelization
#' 
#' Two parallel methods are provided to speed up basic GWR algorithm:
#' 
#' - Multithreading (`omp`)
#' - NVIDIA GPU Computing (`cuda`)
#' 
#' See the vignettes about parallelization to learn more about this topic.
#'
#' @examples
#' data(LondonHP)
#'
#' # Basic usage
#' gwr_lcr(PURCHASE ~ FLOORSZ + UNEMPLOY, LondonHP, 64, TRUE)
#'
#' # Bandwidth Optimization
#' m <- gwr_lcr(PURCHASE ~ FLOORSZ + UNEMPLOY + PROF, LondonHP, 'AIC', TRUE)
#' m
#' 
#' @seealso `browseVignettes("")`
#'
#' @importFrom stats na.action model.frame model.extract model.matrix terms
#' @export
gwr_lcr <- function(
    formula,
    data,
    bw = NA,
    adaptive = FALSE,
    kernel = c("gaussian", "exp", "bisquare", "tricube", "boxcar"),
    longlat = FALSE,
    p = 2.0,
    theta = 0.0,
    lambda = 0,
    lambda_adjust = FALSE,
    cn_thresh = 30,
    cv = TRUE,
    hatmatrix = TRUE,
    parallel_method = c("no", "omp"),
    parallel_arg = c(0),
    verbose = FALSE
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
    mf <- model.frame(formula = formula(formula), data = sf::st_drop_geometry(data))
    mt <- attr(mf, "terms")
    y <- model.extract(mf, "response")
    x <- model.matrix(mt, mf)
    dep_var <- as.character(attr(terms(mf), "variables")[[2]])
    has_intercept <- attr(terms(mf), "intercept") == 1
    indep_vars <- colnames(x)
    indep_vars[which(indep_vars == "(Intercept)")] <- "Intercept"
    colnames(x) <- indep_vars
    if (has_intercept && indep_vars[1] != "Intercept") {
        stop("Please put Intercept to the first column.")
    }

    ### Check whether bandwidth is valid.
    if (missing(bw)) {
        optim_bw <- TRUE
        bw <- Inf
    } else if (is.numeric(bw) || is.integer(bw)) {
        optim_bw <- FALSE
    } else {
        optim_bw <- TRUE
        if(is.character(bw))
            if(!match.arg(bw, c("CV"))){
                stop("Parameter bw must be 'CV' or a number")
            }
        bw <- Inf
    }

    ### Call solver
    c_result <- tryCatch(gwr_lcr_fit(
        x, y, coords,
        bw, adaptive, enum(kernel), longlat, p, theta,
        has_intercept, hatmatrix,
        enum_list(parallel_method, parallel_types), parallel_arg,
        optim_bw
    ), error = function (e) {
        stop("Error:", conditionMessage(e))
    })
    if (optim_bw)
        bw <- c_result$bandwidth
    betas <- c_result$betas
    yhat <- c_result$yhat
    diagnostic <- c_result$diagnostic
    resi <- y - yhat
    # n_dp <- nrow(coords)
    # rss_gw <- sum(resi * resi)
    # sigma <- rss_gw / (n_dp - 2 * shat_trace[1] + shat_trace[2])
    # betas_se <- sqrt(sigma * betas_se)
    # betas_tv <- betas / betas_se

    ### Create result Layer
    colnames(betas) <- indep_vars
    # colnames(betas_se) <- paste(indep_vars, "SE", sep = ".")
    # colnames(betas_tv) <- paste(indep_vars, "TV", sep = ".")
    sdf_data <- as.data.frame(cbind(
        betas,
        "yhat" = yhat,
        "residual" = resi
        # betas_se,
        # betas_tv
    ))
    sdf_data$geometry <- sf::st_geometry(data)
    sdf <- sf::st_sf(sdf_data)

    ### Return result
    gwlcrm <- list(
        SDF = sdf,
        diagnostic = diagnostic,
        args = list(
            # x = x,
            # y = y,
            # coords = coords,
            bw = bw,
            adaptive = adaptive,
            kernel = kernel,
            longlat = longlat,
            p = p,
            theta = theta,
            hatmatrix = hatmatrix,
            has_intercept = has_intercept,
            parallel_method = parallel_method,
            parallel_arg = parallel_arg,
            optim_bw = optim_bw
        ),
        call = mc,
        indep_vars = indep_vars,
        dep_var = dep_var
    )
    class(gwlcrm) <- "gwlcrm"
    gwlcrm
}



#' Print description of a `gwlcrm` object
#'
#' @param x An `hgwlcrm` object returned by [gwr_lcr()].
#' @param decimal_fmt The format string passing to [base::sprintf()].
#' @inheritDotParams print_table_md
#' 
#' @method print gwlcrm
#' @importFrom stats coef fivenum
#' @rdname print
#' @export
print.gwlcrm <- function(x, decimal_fmt = "%.3f", ...) {
    if (!inherits(x, "gwlcrm")) {
        stop("It's not a gwlcrm object.")
    }

    ### Basic Information
    cat("Results of Ridge Geographically Weighted Regression", fill = T)
    cat("===================================================", fill = T)
    cat("  Formula:", deparse(x$call$formula), fill = T)
    cat("     Data:", deparse(x$call$data), fill = T)
    cat("   Kernel:", x$args$kernel, fill = T)
    cat("Bandwidth:", x$args$bw,
        ifelse(x$args$adaptive, "(Nearest Neighbours)", "(Meters)"),
        ifelse(x$args$optim_bw, paste0(
            "(Optimized accroding to CV)"
        ), ""), fill = T)
    cat("\n", fill = T)

    cat("Summary of Coefficient Estimates", fill = T)
    cat("--------------------------------", fill = T)
    betas <- coef(x)
    beta_fivenum <- t(apply(betas, 2, fivenum))
    colnames(beta_fivenum) <- c("Min.", "1st Qu.", "Median", "3rd Qu.", "Max.")
    rownames(beta_fivenum) <- colnames(betas)
    beta_str <- rbind(
        c("Coefficient", colnames(beta_fivenum)),
        cbind(rownames(beta_fivenum), matrix2char(beta_fivenum, decimal_fmt))
    )
    print_table_md(beta_str, ...)
    cat("\n", fill = T)

    cat("Diagnostic Information", fill = T)
    cat("----------------------", fill = T)
    cat("  RSS:", x$diagnostic$RSS, fill = T)
    cat("  ENP:", x$diagnostic$ENP, fill = T)
    cat("  EDF:", x$diagnostic$EDF, fill = T)
    cat("   R2:", x$diagnostic$RSquare, fill = T)
    cat("R2adj:", x$diagnostic$RSquareAdjust, fill = T)
    cat("  AIC:", x$diagnostic$AIC, fill = T)
    cat(" AICc:", x$diagnostic$AICc, fill = T)
    cat("\n", fill = T)
}

#' @describeIn gwr_lcr Plot the result of basic GWR model.
#'
#' @param x A "gwlcrm" object.
#' @param y Ignored.
#' @param columns Column names to plot.
#'  If it is missing or non-character value, all coefficient columns are plottd.
#' @param \dots Additional arguments passing to [sf::plot()].
#' @method plot gwlcrm
#'
#' @examples
#' plot(m)
#'
#' @export
plot.gwlcrm <- function(x, y, ..., columns) {
    if (!inherits(x, "gwlcrm")) {
        stop("It's not a gwlcrm object.")
    }

    sdf <- sf::st_as_sf(x$SDF)
    sdf_colnames <- names(sf::st_drop_geometry(x$SDF))
    if (!missing(columns) && is.character(columns)) {
        valid_columns <- intersect(columns, sdf_colnames)
        if (length(valid_columns) > 0) {
            sdf <- sdf[valid_columns]
        }
    } else { ### Select coefficient columns.
        sdf <- sdf[x$indep_vars]
    }
    plot(sdf, ...)
}

#' @describeIn gwr_lcr Get coefficients of a basic GWR model.
#'
#' @param object A "gwlcrm" object.
#' @param \dots Additional arguments passing to [coef()].
#' 
#' @examples
#' coef(m)
#' 
#' @method coef gwlcrm
#' @export
coef.gwlcrm <- function(object, ...) {
    if (!inherits(object, "gwlcrm")) {
        stop("It's not a gwlcrm object.")
    }
    sf::st_drop_geometry(object$SDF[object$indep_vars])
}

#' @describeIn gwr_lcr Get fitted values of a basic GWR model.
#'
#' @param object A "gwlcrm" object.
#' @param \dots Additional arguments passing to [fitted()].
#' 
#' @examples
#' fitted(m)
#' 
#' @method fitted gwlcrm
#' @export
fitted.gwlcrm <- function(object, ...) {
    if (!inherits(object, "gwlcrm")) {
        stop("It's not a gwlcrm object.")
    }
    object$SDF[["yhat"]]
}

#' @describeIn gwr_lcr Get residuals of a basic GWR model.
#'
#' @param object A "gwlcrm" object.
#' @param \dots Additional arguments passing to [residuals()].
#' 
#' @examples
#' residuals(m)
#' 
#' @method residuals gwlcrm
#' @export
residuals.gwlcrm <- function(object, ...) {
    if (!inherits(object, "gwlcrm")) {
        stop("It's not a gwlcrm object.")
    }
    object$SDF[["residual"]]
}
