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
    parallel_method = c("none", "omp", "cuda", "cluster"),
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
    hasIntercept <- attr(terms(mf), "intercept") == 1
    indep_vars <- colnames(x)
    indep_vars[which(indep_vars == "(Intercept)")] <- "Intercept"

    ### Check whether bandwidth is valid.
    if (missing(bw))
    {
        optim_bw <- TRUE
        optim_bw_criterion <- "AIC"
        bw <- Inf
    }
    else if (is.numeric(bw) || is.integer(bw)) {
        optim_bw <- FALSE
        optim_bw_criterion <- "AIC"
    }
    else
    {
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
    {
        bw <- c_result$bandwidth
    }
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
    geometry <- sf::st_geometry(data)
    sdf_data$geometry <- geometry
    sdf <- sf::st_sf(sdf_data)

    ### Return result
    gwrm <- list(
        SDF = sdf,
        diagnostic = diagnostic,
        betas = betas,
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
            parallel_method = parallel_method,
            parallel_arg = parallel_arg,
            optim_bw = optim_bw,
            optim_bw_criterion = optim_bw_criterion
        ),
        call = mc
    )
    class(gwrm) <- "gwrm"
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
        stop("It's not a hgwrm object.")
    }

    ### Basic Information
    cat("Geographically Weighted Regression Model", fill = T)
    cat("========================================", fill = T)
    cat("  Formula:", deparse(x$call[[2]]), fill = T)
    cat("     Data:", deparse(x$call[[3]]), fill = T)
    cat("   Kernel:", x$args$kernel, fill = T)
    cat("Bandwidth:", x$args$bw,
        ifelse(x$args$adaptive, "(Nearest Neighbours)", "(Meters)"),
        ifelse(x$args$optim_bw, paste0("(Optimized accroding to ", x$args$optim_bw_criterion, ")"), ""),
        fill = T)
    cat("\n", fill = T)

    cat("Summary of Coefficient Estimates", fill = T)
    cat("--------------------------------", fill = T)
    beta_fivenum <- t(apply(x$betas, 2, fivenum))
    colnames(beta_fivenum) <- c("Min.", "1st Qu.", "Median", "3rd Qu.", "Max.")
    rownames(beta_fivenum) <- colnames(x$betas)
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
