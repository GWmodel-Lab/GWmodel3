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
#' @param hatmatrix If TRUE, great circle will be caculated
#' @param parallel_method Parallel method.
#' @param parallel_arg Parallel method argument.
#'
#' @return A `gwrm` object.
#'
#' @examples
#' data(LondonHP)
#' 
#' gwr_basic(PURCHASE~FLOORSZ+UNEMPLOY, LondonHP, 64, TRUE) # Basic usage
#' gwr_basic(PURCHASE~FLOORSZ+UNEMPLOY, LondonHP, 'AIC', TRUE) # Bandwidth Optimization
#'
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
    mc <- match.call(expand.dots = FALSE)
    mt <- match(c("formula", "data"), names(mc), 0L)
    mf <- mc[c(1L, mt)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    y <- model.extract(mf, "response")
    x <- model.matrix(mt, mf)
    dep_var <- as.character(attr(terms(formula(formula)), which = "variables")[[2]])
    hasIntercept <- attr(terms(mf), "intercept") == 1
    indep_vars <- colnames(x)
    indep_vars[which(indep_vars == "(Intercept)")] <- "Intercept"
    if (hasIntercept && indep_vars[1] != "Intercept") {
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
        optim_bw_criterion <- ifelse(is.character(bw), match.arg(bw, c("CV", "AIC")), "AIC")
        bw <- Inf
    }

    ### Call solver
    c_result <- .c_gwr_basic(
        x, y, coords, bw, adaptive, kernel, longlat, p, theta,
        hatmatrix, hasIntercept, parallel_method, parallel_arg,
        optim_bw, optim_bw_criterion
    )
    if (optim_bw)
        bw <- c_result$bandwidth
    betas <- c_result$betas
    betasSE <- c_result$betasSE
    sTrace <- c_result$sTrace
    sHat <- c_result$sHat
    fitted <- c_result$fitted
    diagnostic <- c_result$diagnostic
    resi <- y - fitted
    n_dp <- nrow(coords)
    rss_gw <- sum(resi * resi)
    sigma <- rss_gw / (n_dp - 2 * sTrace[1] + sTrace[2])
    betasSE <- sqrt(sigma * betasSE)
    betasTV <- betas / betasSE

    ### Create result Layer
    colnames(betas) <- indep_vars
    colnames(betasSE) <- paste(indep_vars, "SE", sep = ".")
    colnames(betasTV) <- paste(indep_vars, "TV", sep = ".")
    sdf_data <- as.data.frame(cbind(betas, "yhat" = fitted, "residual" = resi, betasSE, betasTV))
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
            hasIntercept = hasIntercept,
            parallel_method = parallel_method,
            parallel_arg = parallel_arg,
            optim_bw = optim_bw,
            optim_bw_criterion = optim_bw_criterion
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
#' @param gwrm A "`gwrm`" class object.
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
#'
#' @return A `gwrm` object.
#'
#' @examples
#' data(LondonHP)
#' model_sel(gwr_basic(PURCHASE~FLOORSZ+UNEMPLOY+PROF, LondonHP, 'AIC', TRUE), threshold = 100.0)
#' 
#' @export 
model_sel.gwrm <- function(
    gwrm,
    criterion = c("AIC"),
    threshold = 3.0,
    bw = NA,
    optim_bw = c("no", "AIC", "CV"),
    ...
) {
    if (!inherits(gwrm, "gwrm")) {
        stop("It's not a gwrm object.")
    }
    criterion <- match.arg(criterion)

    ### Check whether bandwidth is valid.
    bw_value <- Inf
    optim_bw <- FALSE
    optim_bw_criterion <- "AIC"
    if (missing(bw) || is.na(bw)) {
        bw_value <- gwrm$args$bw
    } else if (is.numeric(bw) || is.integer(bw)) {
        bw_value <- bw
    } else {
        optim_bw <- TRUE
        optim_bw_criterion <- ifelse(is.character(bw), match.arg(bw, c("CV", "AIC")), "AIC")
    }

    ### Calibrate GWR
    c_result <- with(gwrm$args, .c_gwr_basic(
        x, y, coords, bw_value, adaptive, kernel, longlat, p, theta,
        hatmatrix, hasIntercept, parallel_method, parallel_arg,
        optim_bw, optim_bw_criterion,
        select_model = TRUE, select_model_criterion = criterion,
        select_model_threshold = threshold
    ))
    if (optim_bw) 
        bw <- c_result$bandwidth
    betas <- c_result$betas
    betasSE <- c_result$betasSE
    sTrace <- c_result$sTrace
    sHat <- c_result$sHat
    fitted <- c_result$fitted
    diagnostic <- c_result$diagnostic
    resi <- gwrm$args$y - fitted
    n_dp <- nrow(gwrm$args$coords)
    rss_gw <- sum(resi * resi)
    sigma <- rss_gw / (n_dp - 2 * sTrace[1] + sTrace[2])
    betasSE <- sqrt(sigma * betasSE)
    betasTV <- betas / betasSE

    ### Check select variable names
    indep_vars_selected <- c_result$variables + 1
    if (gwrm$args$hasIntercept)
        indep_vars_selected <- c(1, indep_vars_selected)
    indep_vars <- gwrm$indep_vars[indep_vars_selected]

    ### Create result Layer
    colnames(betas) <- indep_vars
    colnames(betasSE) <- paste(indep_vars, "SE", sep = ".")
    colnames(betasTV) <- paste(indep_vars, "TV", sep = ".")
    sdf_data <- as.data.frame(cbind(betas, "yhat" = fitted, "residual" = resi, betasSE, betasTV))
    sdf_data$geometry <- sf::st_geometry(gwrm$SDF)
    sdf <- sf::st_sf(sdf_data)

    ### Convert model sel criterions
    model_sel_criterions <- c_result$model_sel_criterions
    model_sel_criterions$models <- lapply(c_result$model_sel_criterions$models, function(model_var_idx, indep_vars, dep_var, has_intercept) {
        sel_indep_vars <- indep_vars[model_var_idx + 1]
        if (!has_intercept)
            sel_indep_vars <- c("0", sel_indep_vars)
        paste(dep_var, paste(sel_indep_vars, collapse = "+"), sep = "~")
    }, gwrm$indep_vars, gwrm$dep_var, gwrm$args$hasIntercept)
    model_sel_criterions$indep_vars <- gwrm$indep_vars[gwrm$indep_vars != "Intercept"]
    model_sel_criterions$dep_var <- gwrm$dep_var
    class(model_sel_criterions) <- "modelselcritl"

    ### Return result
    gwrm$SDF <- sdf
    gwrm$args$bw <- bw_value
    gwrm$args$select_model <- TRUE
    gwrm$args$select_model_criterion <- criterion
    gwrm$args$select_model_threshold <- threshold
    gwrm$model_sel <- model_sel_criterions
    gwrm$indep_vars <- indep_vars
    gwrm
}

#' Print description of a `gwrm` object
#'
#' @param x An `hgwrm` object returned by [gwr_basic()].
#' @param decimal.fmt The format string passing to [base::sprintf()].
#' @inheritDotParams print_table_md
#'
#' @return No return.
#'
#' @export
#'
print.gwrm <- function(x, decimal.fmt = "%.3f", ...) {
    if (!inherits(x, "gwrm")) {
        stop("It's not a gwrm object.")
    }

    ### Basic Information
    cat("Geographically Weighted Regression Model", fill = T)
    cat("========================================", fill = T)
    cat("  Formula:", deparse(formula(paste(x$dep_var, paste(x$indep_vars, collapse = "+"), sep = "~"))), fill = T)
    cat("     Data:", deparse(x$call[[3]]), fill = T)
    cat("   Kernel:", x$args$kernel, fill = T)
    cat("Bandwidth:", x$args$bw,
        ifelse(x$args$adaptive, "(Nearest Neighbours)", "(Meters)"),
        ifelse(x$args$optim_bw, paste0("(Optimized accroding to ", x$args$optim_bw_criterion, ")"), ""),
        fill = T)
    cat("\n", fill = T)

    cat("Summary of Coefficient Estimates", fill = T)
    cat("--------------------------------", fill = T)
    betas <- sf::st_drop_geometry(x$SDF[, x$indep_vars])
    beta_fivenum <- t(apply(betas, 2, fivenum))
    colnames(beta_fivenum) <- c("Min.", "1st Qu.", "Median", "3rd Qu.", "Max.")
    rownames(beta_fivenum) <- colnames(betas)
    beta_str <- rbind(
        c("Coefficient", colnames(beta_fivenum)),
        cbind(rownames(beta_fivenum), matrix2char(beta_fivenum, decimal.fmt))
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
