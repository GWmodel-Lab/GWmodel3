#' Geographically Weighted Density Regression
#'
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
        config <- rep(config, times = ncol(x))
    } else if (length(config) != ncol(x)) {
        stop("The length of config mush be equal to the number of independent variables.")
    }
    if (has_intercept) {
        config[[1]]@centered <- FALSE
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

    c_result <- .c_gwdr_fit(
        x, y, coords, bw_value, adaptive, kernel,
        has_intercept, TRUE, parallel_method, parallel_arg,
        optim_bw, optim_bw_criterion, optim_bw_threshold,
        optim_bw_step, optim_bw_max_iter
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


