#' Geographically and Temporally Weighted Regression
#'
#' @param formula Formula for regression.
#' @param data A `sf` objects.
#' @param timestamps A vector timestamps for all samples.
#'  Either a numerical vector, or a character vector to be parsed according to `time_format`.
#' @param bw Either a value to set the size of bandwidth,
#'  or one of the following characters to set the criterion for
#'  bandwidth auto-optimization process.
#'  - `AIC`
#'  - `CV`
#'  Note that if `NA` or other non-numeric value is setted,
#'  this parameter will be reset to `Inf`.
#' @param lambda Either a value between 0 and 1 as the weight of temporal distance,
#'  or a `NA` value to enable auto-selection.
#' @param adaptive Whether the bandwidth value is adaptive or not.
#' @param kernel Kernel function used.
#' @param longlat Whether the coordinates
#' @param p Power of the Minkowski distance,
#'  default to 2, i.e., Euclidean distance.
#' @param theta Angle in radian to roate the coordinate system, default to 0.
#' @param optim_bw_range Bounds on bandwidth optimization, a vector of two numeric elements.
#'  Set to `NA_real_` to enable default values selected by the algorithm.
#' @param hatmatrix If TRUE, great circle will be caculated.
#' @param time_format The format used to parse timestamps if they are characters.
#' @param parallel_method Parallel method.
#' @param parallel_arg Parallel method argument.
#' @param verbose Whether to print additional information.
#'
#' @return A `gtwrm` object.
#' 
#' @details
#' ## Parallelization
#' 
#' Supported method(s):
#' 
#' - Multithreading (`omp`)
#' 
#' See the vignettes about parallelization to learn more about this topic.
#'
#' @examples
#' data(LondonHP)
#' LondonHP$time <- as.integer(round(runif(nrow(LondonHP), 1, 365)))
#'
#' # Basic usage
#' m <- gtwr(PURCHASE~FLOORSZ+UNEMPLOY, LondonHP, "time", 64, 0.05, TRUE)
#' m
#' head(coef(m))
#' head(fitted(m))
#' head(residuals(m))
#'
#' # Bandwidth Optimization
#' gtwr(PURCHASE~FLOORSZ+UNEMPLOY, LondonHP, LondonHP$time, NA, 0.05, TRUE)
#' 
#' # Lambda Optimization
#' gtwr(PURCHASE~FLOORSZ+UNEMPLOY, LondonHP, LondonHP$time, 64, NA, TRUE)
#' 
#' # Bandwidth and Lambda optimization
#' gtwr(PURCHASE~FLOORSZ+UNEMPLOY, LondonHP, LondonHP$time, NA, NA, TRUE)
#'
#' @importFrom stats na.action model.frame model.extract model.matrix terms
#' @export
gtwr <- function(
    formula,
    data,
    timestamps,
    bw = NA,
    lambda = NA,
    adaptive = FALSE,
    kernel = c("gaussian", "exp", "bisquare", "tricube", "boxcar"),
    longlat = FALSE,
    p = 2.0,
    theta = 0.0,
    optim_bw_range = c(0, Inf),
    hatmatrix = TRUE,
    time_format = NA,
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

    if (missing(timestamps)) {
        stop("Timestamps are required.")
    } else if (is.character(timestamps)) {
        if (length(timestamps) == 1) {
            if (timestamps %in% names(data)) {
                timestamps <- data[[timestamps]]
            } else {
                stop("The specified `timestamps` must be a (only one) valid column name in `data`.")
            }
        } else if (length(timestamps) == nrow(coords)) {
            timestamps <- as.numeric(as.Date(timestamps, format = time_format))
        } else {
            stop("If `timestamps` is characters, it must be either a column name in `data`, or date-time character vector for all samples.")
        }
    } else if (is.numeric(timestamps) || is.integer(timestamps)) {
        if (length(timestamps) != nrow(coords)) {
            stop("The length of `timestamps` must be equal to the number of samples.")
        }
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

    if (missing(lambda)) {
        optim_lambda <- TRUE
        lambda <- 0.05
    } else if (is.numeric(lambda) && lambda >= 0 && lambda <= 1) {
        optim_lambda <- FALSE
    } else {
        optim_lambda <- TRUE
        lambda <- 0.05
    }

    ### Call solver
    c_result <- tryCatch(gtwr_fit(
        x, y, coords, timestamps, bw, adaptive, enum(kernel), lambda, longlat, p, theta,
        has_intercept, hatmatrix,
        enum_list(parallel_method, parallel_types), parallel_arg,
        optim_bw, enum(optim_bw_criterion, c("AIC", "CV")),
        optim_lambda, as.integer(verbose)
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
    gtwrm <- list(
        SDF = sdf,
        diagnostic = diagnostic,
        args = list(
            x = x,
            y = y,
            coords = coords,
            timestamps = timestamps,
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
            optim_lambda = optim_lambda,
            verbose = verbose
        ),
        call = mc,
        indep_vars = indep_vars,
        dep_var = dep_var
    )
    class(gtwrm) <- c("gtwrm", "gwrm")
    gtwrm
}

#' @describeIn print.gwrm
#'
#' @param x Object returned by GW modelling methods.
#' @param decimal_fmt The format string passing to [base::sprintf()].
#' @inheritDotParams print_table_md
#' 
#' @method print gtwrm
#' @importFrom stats coef fivenum
#' @rdname print
#' @export
print.gtwrm <- function(x, decimal_fmt = "%.3f", ...) {
    if (!inherits(x, "gtwrm")) {
        stop("It's not a gtwrm object.")
    }

    ### Basic Information
    cat("Geographically and Temporally Weighted Regression Model", fill = T)
    cat("=======================================================", fill = T)
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

#' @describeIn gtwr [Not implemented] Predict on new locations.
#'
#' @param object A "gtwrm" object.
#' @param regression_points Data of new locations.
#' @param \dots Additional arguments.
#' @param verbose Whether to print additional message.
#' 
#' @method predict gtwrm
#' @noRd
#' @export
predict.gtwrm <- function(object, regression_points, verbose = FALSE, ...) {
    if (!inherits(object, "gtwrm")) {
        stop("It's not a gtwrm object.")
    }
    stop("GTWR is not predictable now.")
}
