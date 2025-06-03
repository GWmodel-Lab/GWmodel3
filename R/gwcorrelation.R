#' Geographically Weighted Correlation
#' 
#' @param formula A formula, seperated by '~' and joined by '+'
#' @param data A `sf` objects.
#' @param config Parameter-specified weighting configuration.
#'  It must be a list of `GWCorr` objects.
#'  Please find more details in the details section.
#' @param parallel_method Parallel method.
#' @param parallel_arg Parallel method argument.
#' 
#' @return A `gwcorrm` object.
#' 
#' @examples
#' data(LondonHP)
#' # Basic use
#' m <- gwcorrelation(formula = PURCHASE ~ FLOORSZ + UNEMPLOY, data = LondonHP)
#' # Specify more variable pairs
#' gwcorrelation(
#'  formula = PURCHASE + FLOORSZ ~ UNEMPLOY + PROF,
#'  data = LondonHP
#' )
#' # Specify using default for all variables
#' gwcorrelation(
#'   formula = PURCHASE + FLOORSZ ~ UNEMPLOY + PROF,
#'   data = LondonHP, 
#'   config = list(.default = gwcorr_config(adaptive = TRUE, kernel = "gaussian", optim_bw="CV"))
#' )
#' # or
#' gwcorrelation(
#'  formula = PURCHASE + FLOORSZ ~ UNEMPLOY + PROF,
#'  data = LondonHP,
#'  config = list(gwcorr_config(adaptive = FALSE, kernel = "bisquare"))
#' )
#'
#' # Specify more configurations for each variables
#' # Make sure the numbers of variables are correct
#' gwcorrelation(
#'  formula = PURCHASE + FLOORSZ ~ UNEMPLOY + PROF,
#'  data = LondonHP,
#'  config = list(
#'      gwcorr_config(adaptive = TRUE, kernel = "gaussian"),
#'      gwcorr_config(adaptive = TRUE, kernel = "gaussian"),
#'      gwcorr_config(adaptive = TRUE, kernel = "bisquare"),
#'      gwcorr_config(adaptive = TRUE, kernel = "bisquare")
#'  ))
#' # or names should be paired with '_'
#' gwcorrelation(
#'  formula = PURCHASE + FLOORSZ ~ UNEMPLOY + PROF,
#'  data = LondonHP,
#'  config = list(
#'      PURCHASE_UNEMPLOY = gwcorr_config(kernel="bisquare", adaptive = TRUE),
#'      PURCHASE_PROF = gwcorr_config(kernel="bisquare", adaptive = FALSE),
#'      FLOORSZ_UNEMPLOY = gwcorr_config(kernel="gaussian", adaptive = TRUE),
#'      FLOORSZ_PROF = gwcorr_config(kernel="gaussian", adaptive = FALSE)
#'  ))
#' 
#' # Specify configurations by variable names and default
#' gwcorrelation(
#'  formula = PURCHASE + FLOORSZ ~ UNEMPLOY + PROF,
#'  data = LondonHP,
#'  config = list(
#'      PURCHASE_UNEMPLOY = gwcorr_config(bw = 20, adaptive = TRUE, kernel = "bisquare"),
#'      .default = gwcorr_config(adaptive = TRUE, kernel = "bisquare")
#'  ))
#' 
#' @export 
gwcorrelation <- function(
    formula,
    data,
    config = list(gwcorr_config()),
    parallel_method = c("no", "omp"),
    parallel_arg = c(0),
    verbose = FALSE
) {

    ### Check args
    if (!inherits(config, "list")) {
        stop("Parmeter config requires a list of GWCorrConfig.")
    } else if (!all(sapply(config, inherits, "GWCorrConfig"))) {
        stop("Each element in parmeter config requires to be GWCorrConfig.")
    } else {
        for (i in config) {
            valid <- validObject(i)
            if (valid != TRUE) {
                stop(valid)
            }
        }
    }
    parallel_method <- match.arg(parallel_method)

    attr(data, "na.action") <- getOption("na.action")

    data <- do.call(na.action(data), args = list(data))
    coords <- as.matrix(sf::st_coordinates(sf::st_centroid(data)))
    if (is.null(coords) || nrow(coords) != nrow(data))
        stop("Missing coordinates.")

    # make new formula
    formula_str_left <- formula[[2]]
    formula_str_right <- paste(deparse(formula[[3]]), collapse = " ")
    if(length(formula_str_left)>1){
        formula_left <- paste0("cbind(", paste(all.vars(formula_str_left), collapse = ", "), ")")
        formula_right <- paste(formula_str_right, " - 1")
        formula_new <- as.formula(paste(formula_left, "~", formula_right))
    }else{
        formula_new <- update(formula, . ~ . - 1)
    }
    # print(formula_new)

    ### Extract variables
    mc <- match.call(expand.dots = FALSE)
    mf <- model.frame(formula = formula(formula_new), data = sf::st_drop_geometry(data))
    mt <- attr(mf, "terms")
    x1 <- model.extract(mf, "response")#y
    if (is.null(dim(x1))) {
        x1 <- matrix(x1, ncol = 1)
    }
    x2 <- model.matrix(mt, mf)
    if(length(formula_str_left)<=1){
        vars1 <- all.vars(formula_str_left)
    }else{
        vars1 <- colnames(x1)
    }
    vars2 <- colnames(x2)
    colnames(x1) <- vars1

    x2 <- as.matrix(x2)
    x1 <- as.matrix(x1)

    # vars_grid <- expand.grid(v2 = vars2, v1 = vars1)
    # vars_all <- paste0(vars_grid$v1, "_", vars_grid$v2)
    vars_all <- as.vector(outer(vars1, vars2, FUN = function(v1, v2) paste0(v1, "_", v2)))

    ### Extract config
    if (is.null(names(config))) {
        #### When configs are not specified by names, deal with it in the old way
        ncols = ifelse(is.null(ncol(x1)),ncol(x2),ncol(x1)*ncol(x2))
        if (length(config) == 1) {
            config <- rep(config, times = ncols)
        } else if (length(config) != ncols) {
            stop("The length of config mush be equal to the number of variables' product.")
        }
    } else {
        #### When configs are specified by names, match it with parameters
        if (all(vars_all %in% names(config))) {
            config <- config[vars_all]
        } else if (".default" %in% names(config)) {
            config <- ifelse(vars_all %in% names(config), config[vars_all], config[".default"])
        } else {
            stop("Either provide configs for all variable combinations using '_', or supply a default config named '.default'.")
        }
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
            optim_bw <- FALSE
            optim_bw_criterion <- x@optim_bw
            initial_type <- "Specified"
        } else {
            bw_value <- Inf
            optim_bw <- TRUE
            optim_bw_criterion <- "CV"
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

    adaptive <- sapply(config, function(x) x@adaptive)
    kernel <- sapply(config, function(x) x@kernel)
    longlat <- sapply(config, function(x) x@longlat)
    p <- sapply(config, function(x) x@p)
    theta <- sapply(config, function(x) x@theta)

    c_results <- tryCatch(gw_correlation_cal(
        x1, x2, coords,
        bw_value, adaptive, enum(kernel, kernel_enums),
        longlat, p, theta,
        enum(initial_type,gwcorr_initial_enums),
        enum(optim_bw_criterion, gwcorr_bw_criterion_enums),
        enum_list(parallel_method, parallel_types),
        parallel_arg,
        vars_all,
        as.integer(verbose)
    ), error = function (e) {
        stop("Error:", conditionMessage(e))
    })

    # print(c_results)

    sdf_data <- data.frame()
    ## copy results
    {
        local_cov <- c_results$Cov
        local_corr <- c_results$Corr
        local_scorr <- c_results$SCorr
        colnames(local_cov) <- paste(vars_all, "LCov", sep = "_")
        colnames(local_corr) <- paste(vars_all, "LCorr", sep = "_")
        colnames(local_scorr) <- paste(vars_all, "LSCorr", sep = "_")
        sdf_data <- as.data.frame(cbind(local_cov, local_corr, local_scorr))
    } 
    
    sdf_data$geometry <- sf::st_geometry(data)
    sdf <- sf::st_sf(sdf_data)

    ### Return result
    gwcorrm <- list(
        SDF = sdf,
        args = list(
            x1 = x1,
            x2 = x2,
            coords = coords,
            bw_init = bw_value,
            bw_optim = c_results$bw_value,
            adaptive = adaptive,
            kernel = kernel,
            longlat = longlat,
            p = p,
            theta = theta,
            initial_type = initial_type,
            optim_bw_criterion = optim_bw_criterion,
            parallel_method = parallel_method,
            parallel_arg = parallel_arg
        ),
        call = mc,
        names_x1 = vars1,
        names_x2 = vars2,
        names_all = vars_all
    )
    class(gwcorrm) <- "gwcorrm"
    gwcorrm
}



## old code, no export
gw.correlation <- function(
    data,
    vars1,
    vars2,
    bws = NULL,
    adaptive = FALSE,
    kernel = c("gaussian", "exp", "bisquare", "tricube", "boxcar"),
    approach = c("CV", "AIC"),
    longlat = FALSE,
    p = 2.0,
    theta = 0.0,
    parallel_method = c("no", "omp"),
    parallel_arg = c(0),
    verbose = FALSE
) {
    ### Check args
    kernels = match.arg(kernel)
    parallel_method = match.arg(parallel_method)
    approach = match.arg(approach)
    attr(data, "na.action") <- getOption("na.action")

    ### Extract coords
    data <- do.call(na.action(data), args = list(data))
    coords <- as.matrix(sf::st_coordinates(sf::st_centroid(data)))
    if (is.null(coords) || nrow(coords) != nrow(data))
        stop("Missing coordinates.")

    ### Extract variables
    mc <- match.call(expand.dots = FALSE)
    if (inherits(data, "sf")) {
        df <- sf::st_drop_geometry(data)
    } else {
        stop("Input data should be an sf object.")
    }

    if (!is.character(vars1) || length(vars1) == 0 ) {
        stop("Input 'vars1' should be a character vector")
    }
    if (!is.character(vars2) || length(vars2) == 0 ) {
        stop("Input 'vars2' should be a character vector")
    }

    len.var <- length(vars1) * length(vars2)
    col.nm <- colnames(df)

    var.idx1 <- match(vars1, col.nm)[!is.na(match(vars1, col.nm))]
    if (length(var.idx1) == 0) 
        stop("All variables input doesn't match with data")
    x1 <- df[, var.idx1]
    x1 <- as.matrix(x1)
    colnames(x1) <- vars1

    var.idx2 <- match(vars2, col.nm)[!is.na(match(vars2, col.nm))]
    if (length(var.idx2) == 0) 
        stop("All variables input doesn't match with data")
    x2 <- df[, var.idx2]
    x2 <- as.matrix(x2)
    colnames(x2) <- vars2

    ### Check bandwidth.
    initial_type <- c(rep("Specified", length(bws)), rep("Null", ifelse(len.var>=length(bws),len.var - length(bws),len.var)))
    bw_select_type <- c(rep("Specified", length(bws)), rep("AutoSelect", ifelse(len.var>=length(bws),len.var - length(bws),len.var)))
    if (missing(approach))  approach <- "AIC"
    optim_bw_criterion <- rep(match.arg(approach), len.var)

    if (missing(bws)) {
        bws <- c(rep(Inf, len.var))
    } else {
        if (is.numeric(bws) || is.integer(bws)) {
            if (any(bws <= 0)) stop("Bandwidth must be positive numbers.")
        } else {
            stop("Bandwidth must be a number.")
        }
    }

    if (length(bws) > len.var) {
        bws <- bws[1:len.var]
    }

    vars_all <- as.vector(outer(vars1, vars2, FUN = function(v1, v2) paste0(v1, "_", v2)))

    vkernel <- c(rep(kernel, len.var))
    vadaptive <- c(rep(adaptive, len.var))
    vlonglat <- c(rep(longlat, len.var))
    vp <- c(rep(p, len.var))
    vtheta <- c(rep(theta, len.var))

    c_results <- tryCatch(gw_correlation_cal(
        x1, x2, coords,
        bws, vadaptive, enum(vkernel, kernel_enums),
        vlonglat, vp, vtheta,
        enum(initial_type,gwcorr_initial_enums), enum(optim_bw_criterion,gwcorr_bw_criterion_enums),
        enum_list(parallel_method, parallel_types), parallel_arg,
        vars_all, as.integer(verbose)
    ), error = function (e) {
        stop("Error:", conditionMessage(e))
    })

    sdf_data <- data.frame()
    {
        local_cov <- c_results$Cov
        local_corr <- c_results$Corr
        local_scorr <- c_results$SCorr
        colnames(local_cov) <- paste(vars_all, "LCov", sep = "_")
        colnames(local_corr) <- paste(vars_all, "LCorr", sep = "_")
        colnames(local_scorr) <- paste(vars_all, "LSCorr", sep = "_")
        sdf_data <- as.data.frame(cbind(local_cov, local_corr, local_scorr))
    } 
    sdf_data$geometry <- sf::st_geometry(data)
    sdf <- sf::st_sf(sdf_data)

    gwcorrm <- list(
        SDF = sdf,
        args = list(
            x1 = x1,
            x2 = x2,
            coords = coords,
            bw_init = bws,
            bw_optim = c_results$bw_value,
            adaptive = adaptive,
            kernel = kernel,
            longlat = longlat,
            p = p,
            theta = theta,
            initial_type = bw_select_type,
            optim_bw_criterion = approach,
            parallel_method = parallel_method,
            parallel_arg = parallel_arg
        ),
        call = mc,
        names_x1 = vars1,
        names_x2 = vars2,
        names_all = vars_all
    )
    class(gwcorrm) <- "gwcorrm"
    gwcorrm
}


#' Print description of a `gwcorrm` object
#'
#' @param x An `gwcorrm` object returned by [gwr_basic()].
#' @param decimal_fmt The format string passing to [base::sprintf()].
#' @inheritDotParams print_table_md
#' 
#' @importFrom sf st_drop_geometry
#' @method print gwcorrm
#' @rdname print
#' @export
print.gwcorrm <- function(x, ..., decimal_fmt = "%.3f") {
    if (!inherits(x, "gwcorrm"))
        stop("It's not a 'gwcorrm' object.")
    cat("   ***********************************************************************\n")
    cat("   *                       Package   GWmodel                             *\n")
    cat("   ***********************************************************************\n")
    cat("   ***********************Calibration information*************************\n")
    cat("                   Local summary statistics of GW", fill = T, x$args$mode)
    # cat("   =======================================================================\n", fill = T)
    cat("   Data:", deparse(x$call$data), fill = T)
    cat("   Number of summary points:", nrow(x$args$x1), fill = T)

    cat("\n   Parameter-specified Weighting Configuration", fill = T)
    cat("   -----------------------------------------------------------------------\n", fill = T)

    config_str <- with(x$args, data.frame(
        bw = matrix2char(bw_init, ifelse(adaptive, "%.0f", decimal_fmt)),
        optim_bw = bw_optim,
        unit = ifelse(adaptive, "Nearest Neighbors", "Meters"),
        kernel = kernel,
        longlat = longlat,
        p = p,
        theta = theta,
        criterion = ifelse(bw_optim, optim_bw_criterion, "")
    ))

    config_str$distance <- apply(config_str, 1, function(row) {
        distance_type <- "Euclidean"
        if (row["longlat"] == "TRUE" || row["longlat"] == TRUE) {
            distance_type <- "Geodetic"
        } else if (as.numeric(row["p"]) == 2) {
            distance_type <- "Euclidean"
        } else if (as.numeric(row["p"]) == 1) {
            distance_type <- "Manhattan"
        } else if (is.infinite(as.numeric(row["p"]))) {
            distance_type <- "Chebyshev"
        } else {
            distance_type <- "Generalized Minkowski"
        }
        distance_rotated <- (as.numeric(row["theta"]) != 0 && as.numeric(row["p"]) != 2 && !(row["longlat"] == "TRUE" || row["longlat"] == TRUE))
        if (distance_rotated) {
            paste0(distance_type, " (rotated)")
        } else {
            distance_type
        }
    })

    config_str_print <- config_str[, c("bw", "optim_bw", "unit", "kernel", "criterion", "distance")]

    rownames(config_str_print) <- paste0("   ", x$names_all)
    print.data.frame(config_str_print)
    # cat("\n", fill = T)

    cat("\n   ************************Local Summary Statistics:**********************", fill = T)
    res <- st_drop_geometry(x$SDF)

    df <- as.data.frame(t(sapply(res, summary)))
    output <- capture.output(print(df))
    output <- paste0("   ", output)
    cat(output, sep = "\n")

    # print(as.data.frame(t(sapply(res, summary))))
    cat("   ***********************************************************************\n")
}

#' @describeIn gwss Plot the result of basic GWSS model.
#' 
#' @param x A "gwcorrm" object
#' @param y Ignored
#' @param columnsColumn names to plot.
#'  If it is missing or non-character value, all coefficient columns are plotted.
#' @param \dots Additional arguments passing to [sf::plot()].
#' @method plot gwcorrm
#' 
#' @examples 
#' plot(m)
#' 
#' @export
plot.gwcorrm <- function(x, y, ..., columns) {
    if (!inherits(x, "gwcorrm"))
        stop("It's not a 'gwcorrm' object.")

    sdf <- sf::st_as_sf(x$SDF)
    sdf_colnames <- names(sf::st_drop_geometry(x$SDF))
    if (!missing(columns) && is.character(columns)) {
        valid_columns <- intersect(columns, sdf_colnames)
        if (length(valid_columns) > 0) {
            sdf <- sdf[valid_columns]
        }
    }
    plot(sdf, ...)
}