#' Calibrate a basic GGWR model
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
#' @param hatmatrix If TRUE, great circle will be caculated.
#' @param parallel_method Parallel method.
#' @param parallel_arg Parallel method argument.
#' @param verbose Whether to print additional information.
#'
#' @return A `ggwrm` object.
#' 
#' @details
#' ## Parallelization
#' 
#' One parallel methods are provided to speed up basic GGWR algorithm:
#' 
#' - Multithreading (`omp`)
#' 
#' See the vignettes about parallelization to learn more about this topic.
#'
#' @examples
#' data(LondonHP)
#'
#' # Basic usage
#' gwr_generalized(PURCHASE~FLOORSZ+UNEMPLOY, LondonHP, "poisson", 64, TRUE)
#'
#' # Bandwidth Optimization
#' m <- gwr_generalized(PURCHASE ~ FLOORSZ + UNEMPLOY + PROF, LondonHP, 'AIC', TRUE)
#' m
#' 
#' @seealso `browseVignettes("")`
#'
#' @importFrom stats na.action model.frame model.extract model.matrix terms
#' @export
gwr_generalized <- function(
    formula,
    data,
    family = c("poisson", "binomial"),
    bw = NA,
    adaptive = FALSE,
    kernel = c("gaussian", "exp", "bisquare", "tricube", "boxcar"),
    longlat = FALSE,
    p = 2.0,
    theta = 0.0,
    hatmatrix = TRUE,
    parallel_method = c("no", "omp"),
    parallel_arg = c(0)
) {
    ### Check args
    if(missing(family)){
        family = "poisson"
    }
    family = match.arg(family)
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
        optim_bw_criterion <- "AIC"
        bw <- Inf
    } else if (is.numeric(bw) || is.integer(bw)) {
        optim_bw <- FALSE
        optim_bw_criterion <- "AIC"
    } else {
        optim_bw <- TRUE
        optim_bw_criterion <-
            ifelse(is.character(bw), match.arg(bw, c("CV", "AIC")), "AIC")
        bw <- Inf
    }

    ### Call solver
    c_result <- tryCatch(gwr_generalized_fit(
        x, y, coords, 
        enum(family, c("poisson", "binomial")),
        bw, adaptive, enum(kernel), longlat, p, theta,
        hatmatrix, has_intercept,
        optim_bw, enum(optim_bw_criterion, c("AIC", "CV")),
        enum_list(parallel_method, parallel_types), parallel_arg
    ), error = function (e) {
        stop("Error:", conditionMessage(e))
    })
    if (optim_bw)
        bw <- c_result$bandwidth
    betas <- c_result$betas
    diagnostic <- c_result$diagnostic

    yhat <- c_result$yhat
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
    sdf_data <- data.frame(
        betas,
        yhat = yhat,
        residual = resi
        # betas_se,
        # betas_tv
    )
    sdf_data$geometry <- sf::st_geometry(data)
    sdf <- sf::st_sf(sdf_data)

    ### Return result
    ggwrm <- list(
        SDF = sdf,
        diagnostic = diagnostic,
        args = list(
            x = x,
            y = y,
            coords = coords,
            family = family,
            bw = bw,
            adaptive = adaptive,
            kernel = kernel,
            longlat = longlat,
            p = p,
            theta = theta,
            hatmatrix = hatmatrix,
            has_intercept = has_intercept,
            optim_bw = optim_bw,
            optim_bw_criterion = optim_bw_criterion,
            parallel_method = parallel_method,
            parallel_arg = parallel_arg
        ),
        call = mc,
        indep_vars = indep_vars,
        dep_var = dep_var
    )
    class(ggwrm) <- "ggwrm"
    ggwrm
}

#' Print description of a `ggwrm` object
#'
#' @param x An `ggwrm` object returned by [gwr_generalized()].
#' @param decimal_fmt The format string passing to [base::sprintf()].
#' @inheritDotParams print_table_md
#' 
#' @method print ggwrm
#' @importFrom stats coef fivenum
#' @rdname print
#' @export
print.ggwrm <- function(x, decimal_fmt = "%.3f", ...) {
    if (!inherits(x, "ggwrm")) {
        stop("It's not a ggwrm object.")
    }

    ### Basic Information
    cat("Generalized Geographically Weighted Regression Model", fill = T)
    cat("========================================", fill = T)
    cat("  Formula:", deparse(x$call$formula), fill = T)
    cat("     Data:", deparse(x$call$data), fill = T)
    cat("   Family:", x$args$family, fill = T)
    cat("   Kernel:", x$args$kernel, fill = T)
    cat("Bandwidth:", x$args$bw,
        ifelse(x$args$adaptive, "(Nearest Neighbours)", "(Meters)"),
        ifelse(x$args$optim_bw, paste0(
            "(Optimized accroding to ",
            x$args$optim_bw_criterion,
            ")"
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
    cat("   R2:", x$diagnostic$RSquare, fill = T)
    cat("  AIC:", x$diagnostic$AIC, fill = T)
    cat(" AICc:", x$diagnostic$AICc, fill = T)
    cat("\n", fill = T)
}

#' @describeIn gwr_generalized Plot the result of basic GGWR model.
#'
#' @param x A "ggwrm" object.
#' @param y Ignored.
#' @param columns Column names to plot.
#'  If it is missing or non-character value, all coefficient columns are plottd.
#' @param \dots Additional arguments passing to [sf::plot()].
#' @method plot ggwrm
#'
#' @examples
#' plot(m)
#'
#' @export
plot.ggwrm <- function(x, y, ..., columns) {
    if (!inherits(x, "ggwrm")) {
        stop("It's not a ggwrm object.")
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

#' @describeIn gwr_generalized Get coefficients of a Generalized GWR model.
#'
#' @param object A "ggwrm" object.
#' @param \dots Additional arguments passing to [coef()].
#' 
#' @examples
#' coef(m)
#' 
#' @method coef ggwrm
#' @export
coef.ggwrm <- function(object, ...) {
    if (!inherits(object, "ggwrm")) {
        stop("It's not a ggwrm object.")
    }
    sf::st_drop_geometry(object$SDF[object$indep_vars])
}

# ### 测试问题可能在内核库的predict函数中，应该先设置一下mHasRegressionData，否则导致现在出现矩阵列数计算的错误。
# ### 在内核库中测试通过predict后，再完善这部分内容。
# #' @describeIn gwr_generalized Predict on new locations.
# #'
# #' @param object A "ggwrm" object.
# #' @param regression_points Data of new locations.
# #' @param \dots Additional arguments.
# #' @param verbose Whether to print additional message.
# #'
# #' @method predict ggwrm
# #'
# #' @examples
# #' ### predict(m, LondonHP)
# #'
# #' @export
# predict.ggwrm <- function(object, regression_points, verbose = FALSE, ...) {
#     if (!inherits(object, "ggwrm")) {
#         stop("It's not a ggwrm object.")
#     }
#     if (!inherits(regression_points, c("sf", "sfc", "matrix", "data.frame"))) {
#         stop("The type of regression_points is not supported.")
#     }

#     has_intercept <- object$args$has_intercept
#     ### Extract coordinates (and data, if possible)
#     pcoords <- NULL
#     py <- NULL
#     px <- NULL
#     if (inherits(regression_points, c("sf", "sfc"))) {
#         pcoords <- sf::st_coordinates(sf::st_centroid(regression_points))
#         ### Extract X
#         indep_vars <- object$indep_vars
#         if (has_intercept) indep_vars <- indep_vars[-1]
#         if (all(indep_vars %in% names(regression_points))) {
#             px <- as.matrix(sf::st_drop_geometry(regression_points[indep_vars]))
#         }
#         if (has_intercept) {
#             px <- cbind(1, px)
#         }
#         ### Extract Y
#         dep_var <- object$dep_var
#         if (dep_var %in% names(regression_points)) {
#             py <- regression_points[[dep_var]]
#         }
#     } else if(inherits(regression_points, c("matrix", "data.frame"))) {
#         if (length(dim(regression_points)) == 2 && ncol(regression_points) == 2) {
#             pcoords <- as.matrix(regression_points)
#             if (!is.numeric(pcoords)) {
#                 stop("Cannot convert regression_points to numeric matrix.")
#             }
#         } else {
#             stop("Only matrix or data.frame of 2 columns is supported.")
#         }
#     } else {
#         stop("Cannot extract coordinates of regression_points")
#     }
#     pcoords <- as.matrix(pcoords)

#     ### refer to libgwmodel code GWRGenetalized.cpp:175
#     if(object$args$hatmatrix) warning("hatmatrix is not need in prediction")
#     has_hatmatrix <- FALSE

#     ### Predict coefficients
#     c_betas <- tryCatch(with(object$args, gwr_generalized_predict(
#         pcoords, x, y, coords,
#         enum(family, c("poisson", "binomial")),
#         bw, adaptive, enum(kernel, kernel_enums),
#         longlat, p, theta,
#         has_hatmatrix, has_intercept,
#         enum_list(parallel_method, parallel_types), parallel_arg
#     )), error = function(e) {
#         stop("Error:", conditionMessage(e))
#     })

#     result <- as.data.frame(c_betas)
#     colnames(result) <- object$indep_vars

#     ### If px exists, calculate yhat
#     if (!is.null(px)) {
#         yhat <- rowSums(px * c_betas)
#         result <- cbind(result, "yhat" = yhat)

#         if (!is.null(py)) {
#             resi <- py - yhat
#             result <- cbind(result, "residual" = resi)
#         }
#     }

#     ### If regression_points is a sf object, return sf
#     if (inherits(regression_points, c("sf", "sfc"))) {
#         result$geometry <- sf::st_geometry(regression_points)
#         result <- sf::st_as_sf(result)
#     }

#     result
# }