#' Geographically Weighted Summary Statistics
#'
#' @param formula Regresison model.
#' @param data A `sf` objects.
#' @param bw Bandwidth value
#' @param adaptive Whether the bandwidth value is adaptive or not.
#' @param kernel Kernel function used.
#' @param longlat Whether the coordinates
#' @param p Power of the Minkowski distance,
#'  default to 2, i.e., Euclidean distance.
#' @param theta Angle in radian to roate the coordinate system, default to 0.
#' @param polynomial Polynomial
#' @param hatmatrix Hatmatrix
#'
#' @return A `gwscalem` object.
#'
#' @export
gwr_scalable <- function(
    formula,
    data,
    bw,
    adaptive = FALSE,
    kernel = c("gaussian", "exp", "bisquare", "tricube", "boxcar"),
    longlat = FALSE,
    p = 2.0,
    theta = 0.0,
    polynomial = 3,
    hatmatrix = TRUE)
{
    ### Check args
    kernel <- match.arg(kernel)
    attr(data, "na.action") <- getOption("na.action")

    ### Extract coords
    data <- do.call(na.action(data), args = list(data))
    coords <- as.matrix(sf::st_coordinates(sf::st_centroid(data)))
    if (is.null(coords) || nrow(coords) != nrow(data)) {
        stop("Missing coordinates.")
    }

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

    indep_vars <- colnames(x)

    ## Check whether bandwidth is valid.
    if (missing(bw)) {
        print("missing bw")
        optim_bw <- TRUE
        optim_bw_criterion <- "AIC"
        bw <- Inf
    } else if (is.numeric(bw) || is.integer(bw)) {
        print("numeric bw")
        optim_bw <- FALSE
        optim_bw_criterion <- "AIC"
    } else {
        print("else bw")
        optim_bw <- TRUE
        optim_bw_criterion <- 
            ifelse(is.character(bw), match.arg(bw, c("CV", "AIC")), "AIC")
        print(optim_bw)
        print(optim_bw_criterion)
        bw <- 20
    }


    c_results <- tryCatch(gwr_scalable_fit(
        x, y, coords, bw,
        adaptive, enum(kernel),
        longlat, p, theta,
        polynomial, hatmatrix, has_intercept,
        optim_bw, enum(optim_bw_criterion, c("AIC", "CV"))
    ), error = function(e) {
        stop("Error:", conditionMessage(e))
    })

    sdf_data <- data.frame()

    cv <- c_results$cv
    scale <- c_results$scale
    penalty <- c_results$penalty

    betas <- c_results$betas
    colnames(betas) <- indep_vars

    sdf_data <- as.data.frame(betas)

    sdf_data$geometry <- sf::st_geometry(data)
    sdf <- sf::st_sf(sdf_data)

    ### Return result
    gwscalem <- list(
        SDF = sdf,
        args = list(
            # x = x,
            # coords = coords,
            bw = bw,
            adaptive = adaptive,
            kernel = kernel,
            longlat = longlat,
            p = p,
            theta = theta,
            polynomial = polynomial,
            cv = cv,
            scale = scale,
            penalty = penalty,
            hatmatrix = hatmatrix,
            optim_bw = optim_bw,
            optim_bw_criterion = optim_bw_criterion
        ),
        call = mc,
        indep_vars = indep_vars,
        dep_var = dep_var
    )
    class(gwscalem) <- "gwscalem"
    gwscalem
}

#' Print description of a `gwscalem` object
#'
#' @param x An `gwscalem` object returned by [gwr_basic()].
#' @param decimal_fmt The format string passing to [base::sprintf()].
#' @inheritDotParams print_table_md
#'
#' @method print gwscalem
#' @rdname print
#' @export
print.gwscalem <- function(x, ..., decimal_fmt) {
    if (!inherits(x, "gwscalem")) {
        stop("It's not a 'gwscalem' object.")
    }

    cat("   ***********************************************************************\n")
    cat("   *                         Package   GWmodel3                          *\n")
    cat("   ***********************************************************************\n")
    cat("   *              Results of Geographically Weighted Average             *\n")
    cat("   ***********************************************************************\n")
    cat("\n   *********************Model Calibration Information*********************\n")

    cat("   Formula:", deparse(x$call$formula), fill = T)
    cat("   Data:", deparse(x$call$data), fill = T)
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

    cat("\n   ***********************Local Summary Statistics************************", fill = T)
    res <- st_drop_geometry(x$SDF)
    df <- as.data.frame(t(sapply(res, summary)))
    output <- capture.output(print(df))
    output <- paste0("   ", output)
    cat(output, sep = "\n")
    cat("   ***********************************************************************\n")
}