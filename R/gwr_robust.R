#' Calibrate a Robust GWR model
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
#' @param filter Whether to use filtered algorithm.
#' @param parallel_method Parallel method.
#' @param parallel_arg Parallel method argument.
#' @param verbose Whether to print additional information.
#'
#' @return A `rgwrm` object.
#' 
#' @details
#' ## Parallelization
#' 
#' Multithreading (`omp`) are provided to speed up robust GWR algorithm:
#' 
#' See the vignettes about parallelization to learn more about this topic.
#'
#' @examples
#' data(LondonHP)
#'
#' # Basic usage
#' gwr_robust(PURCHASE ~ FLOORSZ + UNEMPLOY, LondonHP, 64, TRUE)
#'
#' # Bandwidth Optimization
#' m <- gwr_robust(PURCHASE ~ FLOORSZ + UNEMPLOY + PROF, LondonHP, 'AIC', TRUE)
#' m
#' 
#' @seealso `browseVignettes("")`
#'
#' @importFrom stats na.action model.frame model.extract model.matrix terms
#' @export
gwr_robust <- function(
    formula,
    data,
    bw = NA,
    adaptive = FALSE,
    kernel = c("gaussian", "exp", "bisquare", "tricube", "boxcar"),
    longlat = FALSE,
    p = 2.0,
    theta = 0.0,
    optim_bw_range = c(0, Inf),
    hatmatrix = TRUE,
    filter = FALSE,
    parallel_method = c("no", "omp", "cuda"),
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
    c_result <- tryCatch(gwr_robust_fit(
        x, y, coords, bw, adaptive, enum(kernel), longlat, p, theta,
        optim_bw_range[1], optim_bw_range[2],
        hatmatrix, has_intercept, filter,
        enum_list(parallel_method, parallel_types), parallel_arg,
        optim_bw, enum(optim_bw_criterion, c("AIC", "CV")),
        select_model = FALSE, select_model_criterion = 0,
        select_model_threshold = 3.0, indep_vars, as.integer(verbose)
    ), error = function (e) {
        stop("Error:", conditionMessage(e))
    })
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
    rgwrm <- list(
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
            optim_bw_lower = optim_bw_range[1],
            optim_bw_upper = optim_bw_range[2],
            hatmatrix = hatmatrix,
            has_intercept = has_intercept,
            is_filtered = filter,
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
    class(rgwrm) <- c("rgwrm", "gwrm")
    rgwrm
}


#' @describeIn gwr_robust Model selection for Robust GWR model
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
#' @return A `rgwrm` object.
#' @method step rgwrm
#'
#' @examples
#' step(m, threshold = 100.0, bw = Inf, optim_bw = "AIC")
#'
#' @importFrom stats formula
#' @export
step.rgwrm <- function(
    object,
    ...,
    criterion = c("AIC"),
    threshold = 3.0,
    bw = NA,
    optim_bw = c("no", "AIC", "CV")
) {
    # step method in gwr
    res <- NextMethod()
    res$is_filtered <- object$args$is_filtered
    class(res) <- c("rgwrm", "gwrm")
    res
}

#' Print description of a `rgwrm` object
#'
#' @param x An `rgwrm` object returned by [gwr_robust()].
#' @param decimal_fmt The format string passing to [base::sprintf()].
#' @inheritDotParams print_table_md
#' 
#' @method print rgwrm
#' @importFrom stats coef fivenum
#' @rdname print
#' @export
print.rgwrm <- function(x, decimal_fmt = "%.3f", ...) {
    if (!inherits(x, "rgwrm")) {
        stop("It's not a rgwrm object.")
    }

    cat("   ***********************************************************************\n")
    cat("   *                         Package   GWmodel3                          *\n")
    cat("   ***********************************************************************\n")
    cat("   *        Results of Roubust Geographically Weighted Regression        *\n")
    cat("   ***********************************************************************\n")
    cat("\n   *********************Model Calibration Information*********************\n")
    cat("   Formula:", deparse(x$call$formula), fill = T)
    cat("   Data:", deparse(x$call$data), fill = T)
    cat("   Kernel:", x$args$kernel, fill = T)
    cat("   Bandwidth:", x$args$bw,
        ifelse(x$args$adaptive, "(Nearest Neighbours)", "(Meters)"),
        ifelse(x$args$optim_bw, paste0(
            "(Optimized accroding to ",
            x$args$optim_bw_criterion,
            ")"
        ), ""), fill = T)
    cat("   Filtered:", x$args$is_filtered, fill = T)

    cat("\n   ****************Local Summary of Coefficient Estimates*****************", fill = T)
    betas <- coef(x)
    beta_fivenum <- t(apply(betas, 2, fivenum))
    colnames(beta_fivenum) <- c("Min.", "1st Qu.", "Median", "3rd Qu.", "Max.")
    rownames(beta_fivenum) <- paste0("   ", colnames(betas))
    beta_str <- rbind(
        c("   Coefficient", colnames(beta_fivenum)),
        cbind(rownames(beta_fivenum), matrix2char(beta_fivenum, decimal_fmt))
    )
    print_table_md(beta_str, ...)

    cat("\n   ************************Diagnostic Information*************************", fill = T)
    cat("   RSS:", x$diagnostic$RSS, fill = T)
    cat("   ENP:", x$diagnostic$ENP, fill = T)
    cat("   EDF:", x$diagnostic$EDF, fill = T)
    cat("   R2:", x$diagnostic$RSquare, fill = T)
    cat("   R2adj:", x$diagnostic$RSquareAdjust, fill = T)
    cat("   AIC:", x$diagnostic$AIC, fill = T)
    cat("   AICc:", x$diagnostic$AICc, fill = T)
    cat("\n", fill = T)
}