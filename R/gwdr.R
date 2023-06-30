#' Geographically Weighted Density Regression (GWDR)
#'
#' @param formula Regresison model.
#' @param data A `sf` objects.
#' @param config Dimension-specified weighting configuration.
#'  It must be a list of `GWDRConfig` objects.
#'  Its length can be 1 or equal to the number of columns of coordinates in `data`.
#'  When its length is 1, elements will be duplicated for each column in coordinates.
#' @param optim_bw Whether optimize bandwidth after selecting models.
#'  Avaliable values are `no`, `AIC`, and `CV`.
#'  If `no` is specified, the bandwidth specified by argument `bw`
#'  is used in calibrating selected models.
#' @param optim_bw_threshold Threshold of bandwidth optimization.
#' @param optim_bw_step Step size in bandwidth optimization.
#' @param optim_bw_max_iter Maximum of iteration in bandwidth optimization.
#' @param parallel_method Parallel method.
#' @param parallel_arg Parallel method argument.
#' 
#' @return A `gwdrm` object.
#'
#' @examples
#' data(LondonHP)
#' gwdr(PURCHASE~FLOORSZ+UNEMPLOY, LondonHP)
#'
#' ### Specific Bandwidth
#' data(LondonHP)
#' gwdr(PURCHASE~FLOORSZ+UNEMPLOY, LondonHP, list(
#'     gwdr_config(0.2, TRUE, "gaussian"),
#'     gwdr_config(0.2, TRUE, "gaussian")
#' ))
#'
#' ### Optim Bandwidth
#' data(LondonHP)
#' gwdr(PURCHASE~FLOORSZ+UNEMPLOY, LondonHP, list(
#'     gwdr_config(0.618, TRUE, "gaussian")
#' ), optim_bw = "AIC")
#'
#' @importFrom stats model.extract model.matrix terms
#' @export
gwdr <- function(
    formula,
    data,
    config = list(gwdr_config()),
    optim_bw = c("no", "AIC", "CV"),
    optim_bw_threshold = 1e-6,
    optim_bw_step = 0.02,
    optim_bw_max_iter = 1e6,
    parallel_method = c("no", "omp"),
    parallel_arg = c(0)
) {
    ### Check args
    if (!inherits(config, "list")) {
        stop("Parameter config requires a list of GWDRConfig.")
    } else if (!all(sapply(config, inherits, "GWDRConfig"))) {
        stop("Each element in parameter config requires to be GWDRConfig.")
    } else {
        for (i in config) {
            valid <- validObject(i)
            if (valid != TRUE) {
                stop(valid)
            }
        }
    }
    parallel_method <- match.arg(parallel_method)
    optim_bw <- match.arg(optim_bw)

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
    if (length(config) == 1) {
        config <- rep(config, times = ncol(coords))
    } else if (length(config) != ncol(coords)) {
        stop("The length of config mush be equal to the number of coordinate columns.")
    }

    ### Process config
    bw_value <- sapply(config, function(x) {
        ifelse(is.numeric(x@bw) || is.integer(x@bw), x@bw, Inf)
    })
    optim_bw_criterion <- "AIC"
    if (optim_bw == "no") {
        optim_bw <- FALSE
    } else {
        optim_bw_criterion <- optim_bw
        optim_bw <- TRUE
    }
    adaptive <- sapply(config, function(x) x@adaptive)
    kernel <- sapply(config, function(x) x@kernel)

    c_result <- gwdr_fit(
        x, y, coords, bw_value, adaptive, enum(kernel, kernel_enums),
        has_intercept, TRUE,
        enum_list(parallel_method, parallel_types), parallel_arg,
        optim_bw, enum(optim_bw_criterion, gwdr_bw_criterion_enums),
        optim_bw_threshold, optim_bw_step, optim_bw_max_iter
    )
    betas <- c_result$betas
    fitted <- c_result$fitted
    diagnostic <- c_result$diagnostic
    resi <- y - fitted
    if (optim_bw)
    {
        bw_value <- c_result$bw_value
    }

    ### Create result Layer
    colnames(betas) <- indep_vars
    sdf_data <- as.data.frame(cbind(
        betas,
        "yhat" = fitted,
        "residual" = resi
    ))
    sdf_data$geometry <- sf::st_geometry(data)
    sdf <- sf::st_sf(sdf_data)

    ### Return result
    gwdrm <- list(
        SDF = sdf,
        diagnostic = diagnostic,
        call = mc,
        indep_vars = indep_vars,
        dep_var = dep_var,
        args = list(
            x = x,
            y = y,
            coords = coords,
            bw_value = bw_value,
            adaptive = adaptive,
            kernel = kernel,
            optim_bw = optim_bw,
            optim_bw_criterion = optim_bw_criterion,
            optim_bw_threshold = optim_bw_threshold,
            optim_bw_step = optim_bw_step,
            optim_bw_max_iter = optim_bw_max_iter,
            has_intercept = has_intercept,
            parallel_method = parallel_method,
            parallel_arg = parallel_arg
        )
    )
    class(gwdrm) <- "gwdrm"
    gwdrm
}
#' Plot the result of basic GWR model.
#'
#' @param x A "gwdrm" object.
#' @param y Ignored.
#' @param columns Column names to plot.
#'  If it is missing or non-character value, all coefficient columns are plottd.
#' @param \dots Additional arguments passing to [sf::plot()].
#' @method plot gwdrm
#' @name plot
#'
#' @examples
#' data(LondonHP)
#' m <- gwr_basic(PURCHASE~FLOORSZ+UNEMPLOY, LondonHP, 64, TRUE)
#' plot(m)
#'
#' @export
plot.gwdrm <- function(x, y, ..., columns) {
    if (!inherits(x, "gwdrm")) {
        stop("It's not a gwdrm object.")
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

#' Plot the result of basic GWR model.
#'
#' @param x A "gwdrm" object.
#' @param y Ignored.
#' @param columns Column names to plot.
#'  If it is missing or non-character value, all coefficient columns are plottd.
#' @param \dots Additional arguments passing to [sf::plot()].
#' @method plot gwdrm
#' @name plot
#'
#' @examples
#' data(LondonHP)
#' m <- gwr_basic(PURCHASE~FLOORSZ+UNEMPLOY, LondonHP, 64, TRUE)
#' plot(m)
#'
#' @export
plot.gwdrm <- function(x, y, ..., columns) {
    if (!inherits(x, "gwdrm")) {
        stop("It's not a gwdrm object.")
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

#' Model selection for GWDR model
#'
#' @param object A "`gwdrm`" class object.
#' @param criterion The model-selection method.
#'  Currently there is only AIC available.
#' @param threshold The threshold of criterion changes. Default to 3.0.
#' @param config Dimension-specified weighting configuration.
#'  If it is omitted, bandwidth settings will be inherited from `object`.
#'  Otherwise, it must be a list of `GWDRConfig` objects.
#'  Its length can be 1 or equal to the number of columns of coordinates in `data`.
#'  When its length is 1, elements will be duplicated for each column in coordinates.
#' @param optim_bw Whether optimize bandwidth after selecting models.
#'  Avaliable values are `no`, `AIC`, and `CV`.
#'  If `no` is specified, the bandwidth specified by argument `bw`
#'  is used in calibrating selected models.
#' @param optim_bw_threshold Threshold of bandwidth optimization.
#' @param optim_bw_step Step size in bandwidth optimization.
#' @param optim_bw_max_iter Maximum of iteration in bandwidth optimization.
#'
#' @importFrom stats formula
#' @export
model_sel.gwdrm <- function(
    object,
    criterion = c("AIC"),
    threshold = 3.0,
    config = list(gwdr_config()),
    optim_bw = c("no", "AIC", "CV"),
    optim_bw_threshold = 1e-6,
    optim_bw_step = 0.02,
    optim_bw_max_iter = 1e6,
    ...
) {
    if (!inherits(object, "gwdrm")) {
        stop("It's not a gwdrm object.")
    }
    criterion <- match.arg(criterion)
    optim_bw <- match.arg(optim_bw)

    ### Check whether bandwidth is valid.
    bw_value <- object$args$bw_value
    if (!missing(config)) {
        bw_value <- sapply(config, function(x) {
            ifelse(is.numeric(x@bw) || is.integer(x@bw), x@bw, Inf)
        })
    }
    optim_bw_criterion <- "AIC"
    if (optim_bw == "no") {
        optim_bw <- FALSE
    } else {
        optim_bw_criterion <- optim_bw
        optim_bw <- TRUE
    }
    adaptive <- sapply(config, function(x) x@adaptive)
    kernel <- sapply(config, function(x) x@kernel)

    ### Calibrate GWR
    c_result <- with(object$args, .c_gwdr_fit(
        x, y, coords, bw_value, adaptive, kernel,
        has_intercept, hatmatrix = TRUE,
        parallel_method, parallel_arg,
        optim_bw, optim_bw_criterion, optim_bw_threshold,
        optim_bw_step, optim_bw_max_iter,
        select_model = TRUE, select_model_threshold = threshold
    ))
    if (optim_bw)
        bw_value <- c_result$bw_value
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

#' Get coefficients of a basic GWR model.
#'
#' @param object A "gwdrm" object.
#' @param \dots Additional arguments passing to [coef()].
#' @method coef gwdrm
#' @name coef
#'
#' @export
coef.gwdrm <- function(object, ...) {
    if (!inherits(object, "gwdrm")) {
        stop("It's not a gwdrm object.")
    }
    sf::st_drop_geometry(object$SDF[object$indep_vars])
}

#' Get fitted values of a basic GWR model.
#'
#' @param object A "gwdrm" object.
#' @param \dots Additional arguments passing to [fitted()].
#' @method fitted gwdrm
#' @name fitted
#' @export
fitted.gwdrm <- function(object, ...) {
    if (!inherits(object, "gwdrm")) {
        stop("It's not a gwdrm object.")
    }
    object$SDF[["yhat"]]
}

#' Get residuals of a basic GWR model.
#'
#' @param object A "gwdrm" object.
#' @param \dots Additional arguments passing to [residuals()].
#' @method residuals gwdrm
#' @name residuals
#' @export
residuals.gwdrm <- function(object, ...) {
    if (!inherits(object, "gwdrm")) {
        stop("It's not a gwdrm object.")
    }
    object$SDF[["residual"]]
}

#' Print description of a `gwdrm` object
#'
#' @param x An `gwdrm` object returned by [gwdr()].
#' @param decimal_fmt The format string passing to [base::sprintf()].
#' @inheritDotParams print_table_md
#' @return No return.
#' @method print gwdrm
#' @name print
#' 
#' @importFrom stats coef fivenum
#' @export
print.gwdrm <- function(x, decimal_fmt = "%.3f", ...) {
    if (!inherits(x, "gwdrm")) {
        stop("It's not a gwdrm object.")
    }

    ### Basic Information
    cat("Geographically Weighted Density Regression Model", fill = T)
    cat("================================================", fill = T)
    cat("  Formula:", deparse(x$call$formula), fill = T)
    cat("     Data:", deparse(x$call$data), fill = T)

    cat("Dimension-specified Weighting Configuration", fill = T)
    cat("-------------------------------------------", fill = T)
    config_str <- with(x$args, data.frame(
        bw = matrix2char(bw_value, decimal_fmt),
        unit = ifelse(adaptive, "%", "Meters"),
        kernel = kernel
    ))
    rownames(config_str) <- colnames(x$args$coords)
    print.data.frame(config_str)
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
