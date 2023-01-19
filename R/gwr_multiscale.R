#' Multiscale GWR
#'
#' @examples
#' data(LondonHP)
#' m <- gwr_multiscale(
#'  formula = PURCHASE ~ FLOORSZ + UNEMPLOY + PROF,
#'  data = LondonHP
#' )
#'
#' # Specify more configurations for all variables
#' m <- gwr_multiscale(
#'  formula = PURCHASE ~ FLOORSZ + UNEMPLOY + PROF,
#'  data = LondonHP,
#'  config = list(mgwr_config(adaptive = TRUE, kernel = "bisquare"))
#' )
#'
#' # Specify more configurations for each variables
#' m <- gwr_multiscale(
#'  formula = PURCHASE ~ FLOORSZ + UNEMPLOY + PROF,
#'  data = LondonHP,
#'  config = list(
#'      mgwr_config(adaptive = TRUE, kernel = "bisquare"),
#'      mgwr_config(adaptive = TRUE, kernel = "bisquare"),
#'      mgwr_config(adaptive = TRUE, kernel = "bisquare"),
#'      mgwr_config(adaptive = TRUE, kernel = "bisquare")
#'  ))
#'
#' @export
gwr_multiscale <- function(
    formula,
    data,
    config = list(mgwr_config()),
    criterion = c("dCVR", "CVR"),
    hatmatrix = T,
    retry_times = 5,
    max_iterations = 2000,
    parallel_method = c("no", "omp"),
    parallel_arg = c(0)
) {
    ### Check args
    if (!inherits(config, "list")) {
        stop("Parmeter config requires a list of MGWRConfig.")
    } else if (!all(sapply(config, inherits, "MGWRConfig"))) {
        stop("Each element in parmeter config requires to be MGWRConfig.")
    } else {
        for (i in config) {
            valid <- validObject(i)
            if (valid != TRUE) {
                stop(valid)
            }
        }
    }
    parallel_method <- match.arg(parallel_method)
    criterion <- match.arg(criterion)

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
        config <- rep(config, times = ncol(x))
    } else if (length(config) != ncol(x)) {
        stop("The length of config mush be equal to the number of independent variables.")
    }
    if (has_intercept) {
        config[[1]]@centered <- FALSE
    }

    ### Process config
    bw_config <- lapply(config, function(x) {
        if (is.na(x@bw)) {
            bw_value <- Inf
            optim_bw <- TRUE
            optim_bw_criterion <- x@optim_bw
            initial_type <- "Null"
        } else if (is.numeric(x@bw) || is.integer(x@bw)) {
            bw_value <- x@bw
            if (x@optim_bw == "no") {
                optim_bw <- FALSE
                optim_bw_criterion <- "AIC"
                initial_type <- "Specified"
            } else {
                optim_bw <- TRUE
                optim_bw_criterion <- x@optim_bw
                initial_type <- "Initial"
            }
        } else {
            bw_value <- Inf
            optim_bw <- TRUE
            optim_bw_criterion <- "AIC"
            initial_type <- "Null"
        }
        list(
            bw_value = bw_value,
            optim_bw = optim_bw,
            optim_bw_criterion = optim_bw_criterion,
            initial_type = initial_type
        )
    })
    bw_value <- sapply(bw_config, function(x) x$bw_value)
    initial_type <- sapply(bw_config, function(x) x$initial_type)
    optim_bw <- sapply(bw_config, function(x) x$optim_bw)
    optim_bw_criterion <- sapply(bw_config, function(x) x$optim_bw_criterion)
    optim_threshold <- sapply(config, function(x) x@optim_threshold)

    adaptive <- sapply(config, function(x) x@adaptive)
    kernel <- sapply(config, function(x) x@kernel)
    longlat <- sapply(config, function(x) x@longlat)
    p <- sapply(config, function(x) x@p)
    theta <- sapply(config, function(x) x@theta)
    centered <- sapply(config, function(x) x@centered)

    c_result <- .c_gwr_multiscale_fit(
        x, y, coords, bw_value, adaptive, kernel, longlat, p, theta,
        optim_bw, optim_bw_criterion, optim_threshold, initial_type, centered,
        criterion, hatmatrix, has_intercept, retry_times,
        max_iterations, parallel_method, parallel_arg
    )
    bw_value <- c_result$bw_value
    betas <- c_result$betas
    fitted <- c_result$fitted
    diagnostic <- c_result$diagnostic
    resi <- y - fitted

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
    gwrmultiscalem <- list(
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
            longlat = longlat,
            p = p,
            theta = theta,
            optim_bw = optim_bw,
            optim_bw_criterion = optim_bw_criterion,
            optim_threshold = optim_threshold,
            initial_type = initial_type,
            centered = centered,
            criterion = criterion,
            hatmatrix = hatmatrix,
            has_intercept = has_intercept,
            retry_times = retry_times,
            parallel_method = parallel_method,
            parallel_arg = parallel_arg
        )
    )
    class(gwrmultiscalem) <- "gwrmultiscalem"
    gwrmultiscalem
}

#' Print description of a `gwrmultiscalem` object
#'
#' @param x An `gwrmultiscalem` object returned by [gwr_multiscale()].
#' @param decimal_fmt The format string passing to [base::sprintf()].
#' @inheritDotParams print_table_md
#' @return No return.
#' @method print gwrmultiscalem
#' @name print
#'
#' @importFrom stats coef fivenum
#' @export
print.gwrmultiscalem <- function(x, decimal_fmt = "%.3f", ...) {
    if (!inherits(x, "gwrmultiscalem")) {
        stop("It's not a gwrmultiscalem object.")
    }

    ### Basic Information
    cat("Multiscale Geographically Weighted Regression Model", fill = T)
    cat("===================================================", fill = T)
    cat("  Formula:", deparse(x$call$formula), fill = T)
    cat("     Data:", deparse(x$call$data), fill = T)
    cat("\n", fill = T)

    cat("Parameter-specified Weighting Configuration", fill = T)
    cat("-------------------------------------------", fill = T)
    config_str <- with(x$args, data.frame(
        bw = matrix2char(bw_value, ifelse(adaptive, "%i", decimal_fmt)),
        unit = ifelse(adaptive, "NN", "Meters"),
        type = initial_type,
        kernel = kernel,
        longlat = longlat,
        p = p,
        theta = theta,
        optim_bw = optim_bw,
        criterion = ifelse(optim_bw, optim_bw_criterion, ""),
        threshold = ifelse(
            optim_bw,
            matrix2char(optim_threshold, "%e"),
            ""
        ),
        centered = centered
    ))
    rownames(config_str) <- x$indep_vars
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

#' Plot the result of basic GWR model.
#'
#' @param x A "gwrmultiscalem" object.
#' @param y Ignored.
#' @param columns Column names to plot.
#'  If it is missing or non-character value, all coefficient columns are plottd.
#' @param \dots Additional arguments passing to [sf::plot()].
#' @method plot gwrmultiscalem
#' @name plot
#'
#' @examples
#' data(LondonHP)
#' m <- gwr_multiscale(
#'  formula = PURCHASE ~ FLOORSZ + UNEMPLOY + PROF,
#'  data = LondonHP,
#'  config = rep(mgwr_config(36, TRUE, "bisquare", optim_bw = "AIC"), 4)
#' )
#' plot(m)
#'
#' @export
plot.gwrmultiscalem <- function(x, y, ..., columns) {
    if (!inherits(x, "gwrmultiscalem")) {
        stop("It's not a gwrmultiscalem object.")
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

#' Get coefficients of a multiscale GWR model.
#'
#' @param object A "gwrmultiscalem" object.
#' @param \dots Additional arguments passing to [coef()].
#' @method coef gwrmultiscalem
#' @name coef
#'
#' @export
coef.gwrmultiscalem <- function(object, ...) {
    if (!inherits(object, "gwrmultiscalem")) {
        stop("It's not a gwrmultiscalem object.")
    }
    sf::st_drop_geometry(object$SDF[object$indep_vars])
}

#' Get fitted values of a basic GWR model.
#'
#' @param object A "gwrmultiscalem" object.
#' @param \dots Additional arguments passing to [fitted()].
#' @method fitted gwrmultiscalem
#' @name fitted
#' @export
fitted.gwrmultiscalem <- function(object, ...) {
    if (!inherits(object, "gwrmultiscalem")) {
        stop("It's not a gwrmultiscalem object.")
    }
    object$SDF[["yhat"]]
}

#' Get residuals of a basic GWR model.
#'
#' @param object A "gwrmultiscalem" object.
#' @param \dots Additional arguments passing to [residuals()].
#' @method residuals gwrmultiscalem
#' @name residuals
#' @export
residuals.gwrmultiscalem <- function(object, ...) {
    if (!inherits(object, "gwrmultiscalem")) {
        stop("It's not a gwrmultiscalem object.")
    }
    object$SDF[["residual"]]
}
