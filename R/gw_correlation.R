#' Geographically Weighted Correlation
#' 
#' @param data A `sf` objects.
#' @param vars1 a vector of variable names to be calculated.
#' @param vars2 a vector of variable names as the responsed.
#' @param bws Bandwidth value list, auto-select if not provided.
#' @param adaptive Whether the bandwidth value is adaptive or not.
#' @param kernel Kernel function used.
#' @param approach Bandwidth selection.
#' @param longlat Whether the coordinates
#' @param p Power of the Minkowski distance,
#'  default to 2, i.e., Euclidean distance.
#' @param theta Angle in radian to roate the coordinate system, default to 0.
#' @param parallel_method Parallel method.
#' @param parallel_arg Parallel method argument.
#' 
#' @return A `gwcorrm` object.
#' 
#' @examples
#' data(LondonHP)
#' gw_correlation(
#'  formula = PURCHASE ~ FLOORSZ + UNEMPLOY + PROF,
#'  data = LondonHP
#' )
#' # Specify more variable pairs
#' gw_correlation(
#'  formula = PURCHASE + FLOORSZ ~ UNEMPLOY + PROF,
#'  data = LondonHP
#' )
#' 
#' # Specify more configurations for all variables
#' gw_correlation(
#'  formula = PURCHASE + FLOORSZ ~ UNEMPLOY + PROF,
#'  data = LondonHP,
#'  config = list(gwcorr_config(adaptive = TRUE, kernel = "bisquare"))
#' )
#' m
#'
#' # Specify more configurations for each variables
#' gw_correlation(
#'  formula = PURCHASE + FLOORSZ ~ UNEMPLOY + PROF,
#'  data = LondonHP,
#'  config = list(
#'      gwcorr_config(adaptive = TRUE, kernel = "gaussian"),
#'      gwcorr_config(adaptive = TRUE, kernel = "gaussian"),
#'      gwcorr_config(adaptive = TRUE, kernel = "bisquare"),
#'      gwcorr_config(adaptive = TRUE, kernel = "bisquare")
#'  ))
#' 
#' # Specify configurations by variable names, names should be paired with '_'
#' gw_correlation(
#'  formula = PURCHASE + FLOORSZ ~ UNEMPLOY + PROF,
#'  data = LondonHP,
#'  config = list(
#'      PURCHASE_UNEMPLOY = gwcorr_config(bw = 20, adaptive = TRUE, kernel = "bisquare"),
#'      .default = gwcorr_config(adaptive = TRUE, kernel = "bisquare")
#'  ))
#' 
#' @export 
gw_correlation <- function(
    formula,
    data,
    config = list(gwcorr_config()),
    parallel_method = c("no", "omp"),
    parallel_arg = c(0)
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
    x1 <- model.extract(mf, "response")
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

    # vars_grid <- expand.grid(v2 = vars2, v1 = vars1)
    # vars_all <- paste0(vars_grid$v1, "_", vars_grid$v2)
    vars_all <- as.vector(outer(vars1, vars2, FUN = function(v1, v2) paste0(v1, "_", v2)))
    # print(vars_all)

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
            # config <- sapply(vars_all, function(name) {
            #     if (name %in% names(config)) config[[name]] else config[[".default"]]
            # })
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

    adaptive <- sapply(config, function(x) x@adaptive)
    kernel <- sapply(config, function(x) x@kernel)
    longlat <- sapply(config, function(x) x@longlat)
    p <- sapply(config, function(x) x@p)
    theta <- sapply(config, function(x) x@theta)

    c_result <- tryCatch(gw_correlation_cal(
        x1, x2, coords,
        bw_value, adaptive, enum(kernel, kernel_enums),
        longlat, p, theta,
        enum(initial_type,gwcorr_initial_enums),
        enum(optim_bw_criterion, gwcorr_bw_criterion_enums),
        enum_list(parallel_method, parallel_types),
        parallel_arg
    ), error = function (e) {
        stop("Error:", conditionMessage(e))
    })

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
            vars1 = x1,
            vars2 = x2,
            coords = coords,
            bw_init = bw_value,
            bw_optim = c_results$bw_value,
            mode = "Correlation",
            adaptive = adaptive,
            kernel = kernels,
            longlat = longlat,
            p = p,
            theta = theta,
            initial_type = initial_type,
            optim_bw_criterion = optim_bw_criterion,
            parallel_method = parallel_method,
            parallel_arg = parallel_arg
        ),
        call = mc,
        vars1 = vars1,
        vars2 = vars2
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
print.gwcorrm <- function(x, ..., decimal_fmt) {
    if (!inherits(x, "gwcorrm"))
        stop("It's not a 'gwcorrm' object.")
    cat("   ***********************************************************************\n")
    cat("   *                       Package   GWmodel                             *\n")
    cat("   ***********************************************************************\n")
    cat("   ***********************Calibration information*************************\n")
    cat("                   Local summary statistics of GW", fill = T, x$args$mode)
    # cat("   =======================================================================\n", fill = T)
    cat("   Data:", deparse(x$call$data), fill = T)
    cat("   Number of summary points:", nrow(x$args$vars1), fill = T)
    cat("   Kernel function:", x$args$kernel, fill = T)
    switch(x$args$mode,
            "Average" = {
                cat("   Bandwidth:", x$args$bw,
                    ifelse(x$args$adaptive, "(Nearest Neighbours)", "(Meters)"),"\n")
            },
            "Correlation" = {
                var_width <- max(nchar(x$name_vars2))
                bw_matrix <- matrix(x$args$bw, nrow = length(x$name_vars1), ncol = length(x$name_vars2), byrow = TRUE)
                cat("   Bandwidth:\n")
                for (i in 1:length(x$name_vars1)) {
                    cat("             ", 
                        sprintf(x$name_vars1[i]),
                        paste(sprintf("%-*s", var_width, x$name_vars2), collapse = " "),fill = TRUE)
                    cat("             ", 
                        paste(sprintf("%-*s", var_width, "")),
                        paste(sprintf("%-*s", var_width, bw_matrix[i, ]), collapse = ""),  
                        ifelse(x$args$adaptive, "(Nearest Neighbours)", "(Meters)"),
                        # ifelse(x$args$optim_bw, paste0("(Optimized according to ", x$args$optim_bw_criterion, ")"), ""),
                        fill = TRUE
                    )
                }
            }
    )
    distance_type <- "Euclidean"
    if (x$args$longlat) distance_type <- "Geodetic"
    else if (x$args$p == 2) distance_type <- "Euclidean"
    else if (x$args$p == 1) distance_type <- "Manhattan"
    else if (is.infinite(x$args$p)) distance_type <- "Chebyshev"
    else distance_type <- "Generalized Minkowski"
    distance_rotated <- (x$args$theta != 0 && x$args$p != 2 && !x$args$longlat)
    cat("   Distance metric:", distance_type, ifelse(distance_rotated, " (rotated)", "distance metric is used."), fill = T)

    cat("\n   ************************Local Summary Statistics:**********************\n", fill = T)
    # cat("   =======================================================================\n", fill = T)
    res <- st_drop_geometry(x$SDF)
    print(as.data.frame(t(sapply(res, summary))))
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