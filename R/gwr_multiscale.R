#' Multiscale GWR
#'
#' @param formula Regresison model.
#' @param data A `sf` objects.
#' @param config Parameter-specified weighting configuration.
#'  It must be a list of `MGWRConfig` objects.
#'  Please find more details in the details section.
#' @param criterion Convergence criterion of back-fitting algorithm.
#' @param optim_bw_range Bounds on bandwidth optimization, a vector of two numeric elements.
#'  Set to `NA_real_` to enable default values selected by the algorithm.
#' @param hatmatrix If TRUE, great circle will be caculated.
#' @param retry_times The number times of continually optimizing each 
#'  parameter-specific bandwidth even though it meets the criterion of convergence,
#'  for avoiding sub-optimal choice due to illusion of convergence.
#' @param max_iterations Maximum number of iterations in the back-fitting procedure.
#' @param parallel_method Parallel method.
#' @param parallel_arg Parallel method argument.
#' @param verbose Output information level.
#'  Can be either `FALSE` or integer values.
#'  Higher values will leads to more output information.
#' 
#' @return A `gwrmultiscalem` object.
#' 
#' @details
#' ## Configuration specification
#' 
#' In the multiscale GWR model,
#' spatial weighting parameters can be specified for each parameter.
#' There are several ways to make it easy and flexible.
#' No matter in which way, the `config` parameter needs to be a list
#' of `MGWRConfig` elements.
#' 
#' When the `config` list is not named, its length needs to be either 1
#' or the number of independent variables (including the intercept if any).
#' For the `config` of length 1,
#' its only value will be applied for every independent variable.
#' For the `config` that as long as independent variables,
#' its values are mapped to variables by position.
#' In other cases, an error will occur to prevent further process.
#' 
#' When the `config` list is named, the names can contain independent-variable names
#' or a special character `".default"`.
#' The function will look up config for each parameter according to its name
#' in the `config` list.
#' If `".default"` can be found in the list,
#' once names of some parameters are missing in the `config`,
#' the function will use the value of name `".default"` instead.
#' However, if not all names can be found in `config` and the `".default"` name is missing,
#' an error will occur to prevent further process.
#' 
#' @examples
#' data(LondonHP)
#' gwr_multiscale(
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
#' m
#'
#' # Specify more configurations for each variables
#' gwr_multiscale(
#'  formula = PURCHASE ~ FLOORSZ + UNEMPLOY + PROF,
#'  data = LondonHP,
#'  config = list(
#'      mgwr_config(adaptive = TRUE, kernel = "bisquare"),
#'      mgwr_config(adaptive = TRUE, kernel = "bisquare"),
#'      mgwr_config(adaptive = TRUE, kernel = "bisquare"),
#'      mgwr_config(adaptive = TRUE, kernel = "bisquare")
#'  ))
#' 
#' # Specify configurations by variable names
#' gwr_multiscale(
#'  formula = PURCHASE ~ FLOORSZ + UNEMPLOY + PROF,
#'  data = LondonHP,
#'  config = list(
#'      FLOORSZ = mgwr_config(bw = 20, adaptive = TRUE, kernel = "bisquare"),
#'      .default = mgwr_config(adaptive = TRUE, kernel = "bisquare")
#'  ))
#'
#' @importFrom methods validObject
#' @importFrom stats na.action model.frame model.extract model.matrix terms
#' @export
gwr_multiscale <- function(
    formula,
    data,
    config = list(mgwr_config()),
    criterion = c("CVR", "dCVR"),
    optim_bw_range = c(0, Inf),
    hatmatrix = T,
    retry_times = 5,
    max_iterations = 2000,
    parallel_method = c("no", "omp"),
    parallel_arg = c(0),
    verbose = FALSE
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

    ### Extract config
    if (is.null(names(config))) {
        #### When configs are not specified by names, deal with it in the old way
        if (length(config) == 1) {
            config <- rep(config, times = ncol(x))
        } else if (length(config) != ncol(x)) {
            stop("The length of config mush be equal to the number of independent variables.")
        }
    } else {
        #### When configs are specified by names, match it with parameters
        if (all(indep_vars %in% names(config))) {
            config <- config[indep_vars]
        } else if (".default" %in% names(config)) {
            config <- ifelse(indep_vars %in% names(config), config[indep_vars], config[".default"])
        } else {
            stop("Either specific configs for all variables, or provide a default config naming '.default'!")
        }
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

    c_result <- tryCatch(gwr_multiscale_fit(
        x, y, coords,
        bw_value, adaptive, enum(kernel, kernel_enums),
        longlat, p, theta,
        optim_bw, enum(optim_bw_criterion, mgwr_bw_criterion_enums),
        optim_threshold, enum(initial_type, mgwr_initial_enums), centered,
        optim_bw_range[1], optim_bw_range[2],
        enum(criterion), hatmatrix, has_intercept, retry_times, max_iterations,
        enum_list(parallel_method, parallel_types), parallel_arg,
        indep_vars, as.integer(verbose)
    ), error = function (e) {
        stop("Error:", conditionMessage(e))
    })
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
            optim_bw_lower = optim_bw_range[1],
            optim_bw_upper = optim_bw_range[2],
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
#'
#' @method print gwrmultiscalem
#' @importFrom stats coef fivenum
#' @rdname print
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
    cat(" AICc:", x$diagnostic$AICc, fill = T)
    cat("\n", fill = T)
}

#' @describeIn gwr_multiscale Plot the result of basic GWR model.
#'
#' @param x A "gwrmultiscalem" object.
#' @param y Ignored.
#' @param columns Column names to plot.
#'  If it is missing or non-character value, all coefficient columns are plottd.
#' @param \dots Additional arguments passing to [sf::plot()].
#' @method plot gwrmultiscalem
#'
#' @examples
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

#' @describeIn gwr_multiscale Get coefficients of a multiscale GWR model.
#'
#' @param object A "gwrmultiscalem" object.
#' @param \dots Additional arguments passing to [coef()].
#'
#' @examples 
#' coef(m)
#'
#' @method coef gwrmultiscalem
#' @export
coef.gwrmultiscalem <- function(object, ...) {
    if (!inherits(object, "gwrmultiscalem")) {
        stop("It's not a gwrmultiscalem object.")
    }
    sf::st_drop_geometry(object$SDF[object$indep_vars])
}

#' @describeIn gwr_multiscale Get fitted values of a basic GWR model.
#'
#' @param object A "gwrmultiscalem" object.
#' @param \dots Additional arguments passing to [fitted()].
#'
#' @examples
#' fitted(m)
#'
#' @method fitted gwrmultiscalem
#' @export
fitted.gwrmultiscalem <- function(object, ...) {
    if (!inherits(object, "gwrmultiscalem")) {
        stop("It's not a gwrmultiscalem object.")
    }
    object$SDF[["yhat"]]
}

#' @describeIn gwr_multiscale Get residuals of a basic GWR model.
#'
#' @param object A "gwrmultiscalem" object.
#' @param \dots Additional arguments passing to [residuals()].
#'
#' @examples 
#' residuals(m)
#'
#' @method residuals gwrmultiscalem
#' @export
residuals.gwrmultiscalem <- function(object, ...) {
    if (!inherits(object, "gwrmultiscalem")) {
        stop("It's not a gwrmultiscalem object.")
    }
    object$SDF[["residual"]]
}
