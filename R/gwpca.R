#' Calibrate a GWPCA model
#'
#' @param formula Regresison model.
#' @param data A `sf` objects.
#' @param bw Should provide a value to set the size of bandwidth,
#' @param adaptive Whether the bandwidth value is adaptive or not.
#' @param kernel Kernel function used.
#' @param longlat Whether the coordinates
#' @param p Power of the Minkowski distance,
#'  default to 2, i.e., Euclidean distance.
#' @param theta Angle in radian to roate the coordinate system, default to 0.
#' @param components How many components want to keep, default to 2.
#'
#' @return A `gwpcam` object.
#'
#' @examples
#' data(LondonHP)
#'
#' # Basic usage
#' gwpca(~FLOORSZ + UNEMPLOY + PROF, LondonHP, 50, TRUE, components = 2)
#'
#' @seealso `browseVignettes("")`
#'
#' @importFrom stats na.action model.frame model.extract model.matrix terms
#' @export
gwpca <- function(
    formula,
    data,
    bw = NA,
    adaptive = FALSE,
    kernel = c("gaussian", "exp", "bisquare", "tricube", "boxcar"),
    longlat = FALSE,
    p = 2.0,
    theta = 0.0,
    components = 2
) {
    ### Check args
    kernel = match.arg(kernel)
    attr(data, "na.action") <- getOption("na.action")

    ### Extract coords
    data <- do.call(na.action(data), args = list(data))
    coords <- as.matrix(sf::st_coordinates(sf::st_centroid(data)))
    if (is.null(coords) || nrow(coords) != nrow(data))
        stop("Missing coordinates.")

    ### Extract variables
    mc <- match.call(expand.dots = FALSE)

    mf <- model.frame(formula = formula(update(formula, ~ . - 1)), data = sf::st_drop_geometry(data))
    mt <- attr(mf, "terms")
    x <- model.matrix(mt, mf)
    indep_vars <- colnames(x)

    ### Check whether bandwidth is valid.
    if (missing(bw)) {
        stop("Please input bandwidth!")
    }
    if (!(is.numeric(bw) || is.integer(bw))) {
        stop("Please input correct bandwidth!")
    }

    components <- as.integer(floor(components))
    if (components < 1) {
        stop("Components must be an interger greater than 0 !!")
    }
    if(length(indep_vars) < components) {
        stop("Components to keep must be lower than variable counts!")
    }



    c_result <- tryCatch(gwpca_cal(
        x, coords,
        bw, adaptive, enum(kernel), longlat, p, theta,
        components
    ), error = function(e) {
        stop("Error:", conditionMessage(e))
    })

    local_loadings <- c_result$local_loadings
    # local_scores <- c_result$local_scores

    local_PV <- c_result$local_PV

    pc_names <- paste0("PC", seq_len(components), "_PV")
    colnames(local_PV) <- pc_names
    sdf_data <- as.data.frame(local_PV)
    sdf_data$geometry <- sf::st_geometry(data)
    sdf <- sf::st_sf(sdf_data)

    ### Return result
    gwpcam <- list(
        SDF = sdf,
        local_loadings = local_loadings,
        # local_scores = local_scores,
        args = list(
            x = x,
            coords = coords,
            bw = bw,
            adaptive = adaptive,
            kernel = kernel,
            longlat = longlat,
            p = p,
            theta = theta,
            keep_components = components
        ),
        call = mc,
        indep_vars = indep_vars
    )
    class(gwpcam) <- "gwpcam"
    gwpcam
}

#' Print description of a `gwpcam` object
#'
#' @param x An `gwpcam` object returned by [gwpca()].
#' @param decimal_fmt The format string passing to [base::sprintf()].
#' @inheritDotParams print_table_md
#'
#' @method print gwpcam
#' @importFrom stats coef fivenum
#' @rdname print
#' @export
print.gwpcam <- function(x, decimal_fmt = "%.3f", ...) {
    if (!inherits(x, "gwpcam")) {
        stop("It's not a gwpcam object.")
    }

    ### Basic Information
    cat("Geographically Weighted Principal Component Analysis", fill = T)
    cat("========================================", fill = T)
    cat("  Formula:", deparse(x$call$formula), fill = T)
    cat("     Data:", deparse(x$call$data), fill = T)
    cat("   Kernel:", x$args$kernel, fill = T)
    cat("Bandwidth:", x$args$bw, ifelse(x$args$adaptive, "(Nearest Neighbours)", "(Meters)"), fill = T)
    cat("   Keep components:", x$args$keep_components, fill = T)
    cat("\n", fill = T)

    cat("Summary of Local PV", fill = T)
    cat("--------------------------------", fill = T)
    local_PV <- sf::st_drop_geometry(x$SDF)
    local_PV_fivenum <- t(apply(local_PV, 2, fivenum))
    colnames(local_PV_fivenum) <- c("Min.", "1st Qu.", "Median", "3rd Qu.", "Max.")
    rownames(local_PV_fivenum) <- colnames(local_PV)
    loadings_1st_str <- rbind(
        c("Local loadings", colnames(local_PV_fivenum)),
        cbind(rownames(local_PV_fivenum), matrix2char(local_PV_fivenum, decimal_fmt))
    )
    print_table_md(loadings_1st_str, ...)
    cat("\n", fill = T)

    cat("Summary of Local Loadings in PC1", fill = T)
    cat("--------------------------------", fill = T)
    loadings_1st <- x$local_loadings[, , 1]
    loadings_1st_fivenum <- t(apply(loadings_1st, 2, fivenum))
    colnames(loadings_1st_fivenum) <- c("Min.", "1st Qu.", "Median", "3rd Qu.", "Max.")
    rownames(loadings_1st_fivenum) <- x$indep_vars
    loadings_1st_str <- rbind(
        c("Local loadings", colnames(loadings_1st_fivenum)),
        cbind(rownames(loadings_1st_fivenum), matrix2char(loadings_1st_fivenum, decimal_fmt))
    )
    print_table_md(loadings_1st_str, ...)
    cat("\n", fill = T)
}