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
#' @param hatmatrix If TRUE, great circle will be caculated.
#' @param parallel_method Parallel method.
#' @param parallel_arg Parallel method argument.
#' @param verbose Whether to print additional information.
#'
#' @return A `gwrm` object.
#'
#' @examples
#' data(LondonHP)
#'
#' # Basic usage
#' gwr_basic(PURCHASE~FLOORSZ+UNEMPLOY, LondonHP, 64, TRUE)
#'
#' # Bandwidth Optimization
#' gwr_basic(PURCHASE~FLOORSZ+UNEMPLOY, LondonHP, 'AIC', TRUE)
#'
#' @importFrom stats model.extract model.matrix terms
#' @export
gwr_basic <- function(
    formula,
    data,
    bw = NA,
    adaptive = FALSE,
    kernel = c("gaussian", "exp", "bisquare", "tricube", "boxcar"),
    longlat = FALSE,
    p = 2.0,
    theta = 0.0,
    hatmatrix = TRUE,
    parallel_method = c("no", "omp"),
    parallel_arg = c(0),
    verbose = FALSE
) {
    ### Check args
    kernel = match.arg(kernel)
    parallel_method = match.arg(parallel_method)

    ### Extract coords
    coords <- as.matrix(sf::st_coordinates(sf::st_centroid(data)))
    if (is.null(coords) || nrow(coords) != nrow(data))
        stop("Missing coordinates.")

    ### Extract variables
    mc <- match.call(expand.dots = FALSE)
    mt <- match(c("formula", "data"), names(mc), 0L)
    mf <- mc[c(1L, mt)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    y <- model.extract(mf, "response")
    x <- model.matrix(mt, mf)
    dep_var <- as.character(attr(terms(formula(formula)), "variables")[[2]])
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
    c_result <- gwr_basic_fit(
        x, y, coords, bw, adaptive, enum(kernel), longlat, p, theta,
        hatmatrix, has_intercept,
        enum_list(parallel_method, parallel_types), parallel_arg,
        optim_bw, enum(optim_bw_criterion, c("AIC", "CV")),
        select_model = FALSE, select_model_criterion = 0,
        select_model_threshold = 3.0, indep_vars, verbose
    )
    if (optim_bw)
        bw <- c_result$bandwidth
    betas <- c_result$betas
    betas_se <- c_result$betasSE
    shat_trace <- c_result$sTrace
    fitted <- c_result$fitted
    diagnostic <- c_result$diagnostic
    resi <- y - fitted
    n_dp <- nrow(coords)
    rss_gw <- sum(resi * resi)
    sigma <- rss_gw / (n_dp - 2 * shat_trace[1] + shat_trace[2])
    betas_se <- sqrt(sigma * betas_se)
    betas_tv <- betas / betas_se

    ### Create result Layer
    colnames(betas) <- indep_vars
    colnames(betas_se) <- paste(indep_vars, "SE", sep = ".")
    colnames(betas_tv) <- paste(indep_vars, "TV", sep = ".")
    sdf_data <- as.data.frame(cbind(
        betas,
        "yhat" = fitted,
        "residual" = resi,
        betas_se,
        betas_tv
    ))
    sdf_data$geometry <- sf::st_geometry(data)
    sdf <- sf::st_sf(sdf_data)

    ### Return result
    gwrm <- list(
        SDF = sdf,
        diagnostic = diagnostic,
        args = list(
            x = x,
            y = y,
            coords = coords,
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
            optim_bw = optim_bw,
            optim_bw_criterion = optim_bw_criterion,
            verbose = verbose
        ),
        call = mc,
        indep_vars = indep_vars,
        dep_var = dep_var
    )
    class(gwrm) <- "gwrm"
    gwrm
}

#' Model selection for basic GWR model
#'
#' @param object A "`gwrm`" class object.
#' @param criterion The model-selection method.
#'  Currently there is only AIC available.
#' @param threshold The threshold of criterion changes. Default to 3.0.
#' @param bw The bandwidth used in selecting models.
#'  If `NA` or missing, the bandwidth value will be derived
#'  from the `gwrm` object.
#' @param optim_bw Whether optimize bandwidth after selecting models.
#'  Avaliable values are `no`, `AIC`, and `CV`.
#'  If `no` is specified, the bandwidth specified by argument `bw`
#'  is used in calibrating selected models.
#' @param \dots Other parameters.
#' @return A `gwrm` object.
#' @method model_sel gwrm
#' @name model_sel
#'
#' @examples
#' data(LondonHP)
#' model_sel(gwr_basic(PURCHASE~FLOORSZ+UNEMPLOY+PROF, LondonHP, 'AIC', TRUE),
#'           threshold = 100.0)
#'
#' @importFrom stats formula
#' @export
model_sel.gwrm <- function(
    object,
    criterion = c("AIC"),
    threshold = 3.0,
    bw = NA,
    optim_bw = c("no", "AIC", "CV"),
    ...
) {
    if (!inherits(object, "gwrm")) {
        stop("It's not a gwrm object.")
    }
    criterion <- match.arg(criterion)

    ### Check whether bandwidth is valid.
    bw_value <- Inf
    if (missing(bw) || is.na(bw)) {
        bw_value <- object$args$bw
        optim_bw_criterion <- "AIC"
        optim_bw <- FALSE
    } else {
        optim_bw_criterion <- match.arg(optim_bw)
        if (optim_bw_criterion == "no") {
            optim_bw_criterion <- "AIC"
            optim_bw <- FALSE
        } else {
            optim_bw <- TRUE
        }
        bw_value <- ifelse(is.numeric(bw) || is.integer(bw), bw, Inf)
    }

    ### Calibrate GWR
    c_result <- with(object$args, gwr_basic_fit(
        x, y, coords, bw_value, adaptive, enum(kernel, kernel_enums),
        longlat, p, theta, hatmatrix, has_intercept,
        enum_list(parallel_method, parallel_types), parallel_arg,
        optim_bw, enum(optim_bw_criterion, c("AIC", "CV")),
        select_model = TRUE, select_model_criterion = enum(criterion),
        select_model_threshold = threshold, object$indep_vars, verbose
    ))
    if (optim_bw)
        bw_value <- c_result$bandwidth
    betas <- c_result$betas
    betas_se <- c_result$betasSE
    shat_trace <- c_result$sTrace
    fitted <- c_result$fitted
    diagnostic <- c_result$diagnostic
    resi <- object$args$y - fitted
    n_dp <- nrow(object$args$coords)
    rss_gw <- sum(resi * resi)
    sigma <- rss_gw / (n_dp - 2 * shat_trace[1] + shat_trace[2])
    betas_se <- sqrt(sigma * betas_se)
    betas_tv <- betas / betas_se

    ### Check select variable names
    indep_vars <- object$indep_vars[c_result$variables + 1]
    if (object$args$has_intercept) {
        formula_up <- reg_formula(object$dep_var, indep_vars)
        indep_vars <- c("Intercept", indep_vars)
    } else {
        formula_up <- reg_formula(object$dep_var, c("0", indep_vars))
    }

    ### Create result Layer
    colnames(betas) <- indep_vars
    colnames(betas_se) <- paste(indep_vars, "SE", sep = ".")
    colnames(betas_tv) <- paste(indep_vars, "TV", sep = ".")
    sdf_data <- as.data.frame(cbind(
        betas,
        "yhat" = fitted,
        "residual" = resi,
        betas_se,
        betas_tv
    ))
    sdf_data$geometry <- sf::st_geometry(object$SDF)
    sdf <- sf::st_sf(sdf_data)

    ### Convert model sel criterions
    model_sel_criterions <- list(
        models = lapply(c_result$model_sel_criterions$models, function(
            model_var_idx,
            indep_vars,
            dep_var,
            has_intercept
        ) {
            sel_indep_vars <- indep_vars[model_var_idx + 1]
            if (!has_intercept)
                sel_indep_vars <- c("0", sel_indep_vars)
            paste(dep_var, paste(sel_indep_vars, collapse = "+"), sep = "~")
        }, object$indep_vars, object$dep_var, object$args$has_intercept),
        criterion_values = c_result$model_sel_criterions$criterions,
        criterion = criterion,
        threshold = threshold,
        indep_vars = object$indep_vars[object$indep_vars != "Intercept"],
        dep_var = object$dep_var
    )
    class(model_sel_criterions) <- "modelselcritl"

    ### Return result
    object$SDF <- sdf
    object$args$x <- object$args$x[, indep_vars]
    object$args$bw <- bw_value
    object$args$select_model <- TRUE
    object$args$select_model_criterion <- criterion
    object$args$select_model_threshold <- threshold
    object$diagnostic <- diagnostic
    object$model_sel <- model_sel_criterions
    object$indep_vars <- indep_vars
    object$call$formula <- str2lang(formula_up)
    object
}

#' Print description of a `gwrm` object
#'
#' @param x An `hgwrm` object returned by [gwr_basic()].
#' @param decimal_fmt The format string passing to [base::sprintf()].
#' @inheritDotParams print_table_md
#' @return No return.
#' @method print gwrm
#' @name print
#' 
#' @importFrom stats coef fivenum
#' @export
print.gwrm <- function(x, decimal_fmt = "%.3f", ...) {
    if (!inherits(x, "gwrm")) {
        stop("It's not a gwrm object.")
    }

    ### Basic Information
    cat("Geographically Weighted Regression Model", fill = T)
    cat("========================================", fill = T)
    cat("  Formula:", deparse(x$call$formula), fill = T)
    cat("     Data:", deparse(x$call$data), fill = T)
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
    cat("  ENP:", x$diagnostic$ENP, fill = T)
    cat("  EDF:", x$diagnostic$EDF, fill = T)
    cat("   R2:", x$diagnostic$RSquare, fill = T)
    cat("R2adj:", x$diagnostic$RSquareAdjust, fill = T)
    cat("  AIC:", x$diagnostic$AIC, fill = T)
    cat(" AICc:", x$diagnostic$AICc, fill = T)
    cat("\n", fill = T)
}

#' Plot the result of basic GWR model.
#'
#' @param x A "gwrm" object.
#' @param y Ignored.
#' @param columns Column names to plot.
#'  If it is missing or non-character value, all coefficient columns are plottd.
#' @param \dots Additional arguments passing to [sf::plot()].
#' @method plot gwrm
#' @name plot
#'
#' @examples
#' data(LondonHP)
#' m <- gwr_basic(PURCHASE~FLOORSZ+UNEMPLOY, LondonHP, 64, TRUE)
#' plot(m)
#'
#' @export
plot.gwrm <- function(x, y, ..., columns) {
    if (!inherits(x, "gwrm")) {
        stop("It's not a gwrm object.")
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

#' Get coefficients of a basic GWR model.
#'
#' @param object A "gwrm" object.
#' @param \dots Additional arguments passing to [coef()].
#' @method coef gwrm
#' @name coef
#'
#' @export
coef.gwrm <- function(object, ...) {
    if (!inherits(object, "gwrm")) {
        stop("It's not a gwrm object.")
    }
    sf::st_drop_geometry(object$SDF[object$indep_vars])
}

#' Get fitted values of a basic GWR model.
#'
#' @param object A "gwrm" object.
#' @param \dots Additional arguments passing to [fitted()].
#' @method fitted gwrm
#' @name fitted
#' @export
fitted.gwrm <- function(object, ...) {
    if (!inherits(object, "gwrm")) {
        stop("It's not a gwrm object.")
    }
    object$SDF[["yhat"]]
}

#' Get residuals of a basic GWR model.
#'
#' @param object A "gwrm" object.
#' @param \dots Additional arguments passing to [residuals()].
#' @method residuals gwrm
#' @name residuals
#' @export
residuals.gwrm <- function(object, ...) {
    if (!inherits(object, "gwrm")) {
        stop("It's not a gwrm object.")
    }
    object$SDF[["residual"]]
}

#' Predict on new locations.
#'
#' @param object A "gwrm" object.
#' @param regression_points Data of new locations.
#' @param \dots Additional arguments.
#' @return A "gwrm" object.
#' @method predict gwrm
#' @name predict
#'
#' @examples
#' data(LondonHP)
#' m <- gwr_basic(PURCHASE ~ FLOORSZ + UNEMPLOY + PROF, LondonHP, 64, TRUE)
#' predict(m, LondonHP)
#'
#' @export
predict.gwrm <- function(object, regression_points, ...) {
    if (!inherits(object, "gwrm")) {
        stop("It's not a gwrm object.")
    }
    if (!inherits(regression_points, c("sf", "sfc", "matrix", "data.frame"))) {
        stop("The type of regression_points is not supported.")
    }

    has_intercept <- object$args$has_intercept
    ### Extract coordinates (and data, if possible)
    pcoords <- NULL
    py <- NULL
    px <- NULL
    if (inherits(regression_points, c("sf", "sfc"))) {
        pcoords <- sf::st_coordinates(sf::st_centroid(regression_points))
        ### Extract X
        indep_vars <- object$indep_vars
        if (has_intercept) indep_vars <- indep_vars[-1]
        if (all(indep_vars %in% names(regression_points))) {
            px <- as.matrix(sf::st_drop_geometry(regression_points[indep_vars]))
        }
        if (has_intercept) {
            px <- cbind(1, px)
        }
        ### Extract Y
        dep_var <- object$dep_var
        if (dep_var %in% names(regression_points)) {
            py <- regression_points[[dep_var]]
        }
    } else if(inherits(regression_points, c("matrix", "data.frame"))) {
        if (length(dim(regression_points)) == 2 && ncol(regression_points) == 2) {
            pcoords <- as.matrix(regression_points)
            if (!is.numeric(pcoords)) {
                stop("Cannot convert regression_points to numeric matrix.")
            }
        } else {
            stop("Only matrix or data.frame of 2 columns is supported.")
        }
    } else {
        stop("Cannot extract coordinates of regression_points")
    }
    pcoords <- as.matrix(pcoords)

    ### Predict coefficients
    c_betas <- with(object$args, gwr_basic_predict(
        pcoords, x, y, coords, bw, adaptive, enum(kernel, kernel_enums),
        longlat, p, theta, has_intercept, 
        enum_list(parallel_method, parallel_types), parallel_arg
    ))

    result <- as.data.frame(c_betas)
    colnames(result) <- object$indep_vars

    ### If px exists, calculate yhat
    if (!is.null(px)) {
        yhat <- rowSums(px * c_betas)
        result <- cbind(result, "yhat" = yhat)

        if (!is.null(py)) {
            resi <- py - yhat
            result <- cbind(result, "residual" = resi)
        }
    }

    ### If regression_points is a sf object, return sf
    if (inherits(regression_points, c("sf", "sfc"))) {
        result$geometry <- sf::st_geometry(regression_points)
        result <- sf::st_as_sf(result)
    }

    result
}
