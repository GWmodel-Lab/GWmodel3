#' GWR Collinearity Diagnostics
#'
#' Compute local collinearity diagnostics including VIF, CN, and VDP
#' for geographically weighted regression. This implementation follows
#' the methodology of GWmodel::gwr.collin.diagno().
#'
#' @param formula Regression model formula.
#' @param data A `sf` object.
#' @param bw Bandwidth value.
#' @param adaptive Whether the bandwidth is adaptive (number of neighbors).
#' @param kernel Kernel function: "gaussian", "exp", "bisquare", "tricube", "boxcar".
#' @param longlat Whether coordinates are longitude/latitude.
#' @param p Power of Minkowski distance (2 = Euclidean).
#' @param theta Angle to rotate coordinates (radians).
#' @param dMat Optional pre-computed distance matrix.
#'
#' @return A list containing:
#' \describe{
#'   \item{SDF}{sf object with VIF, local_CN, VDP, and correlations}
#'   \item{local_CN}{Local condition numbers}
#'   \item{VIF}{Matrix of local VIF values (locations x variables, excluding intercept)}
#'   \item{VDP}{Matrix of variance decomposition proportions for smallest singular value}
#'   \item{corr.mat}{Local correlation coefficients between variable pairs}
#' }
#'
#' @details
#' Collinearity diagnostics interpretation:
#' \itemize{
#'   \item VIF > 10: indicates problematic collinearity
#'   \item CN > 30: indicates severe collinearity
#'   \item VDP: values > 0.5 for two or more variables with high CN indicate collinearity
#' }
#'
#' @references
#' Belsley, D.A., Kuh, E., and Welsch, R.E. (1980). Regression Diagnostics.
#' New York: John Wiley & Sons.
#'
#' Wheeler, D. and Tiefelsdorf, M. (2005). Multicollinearity and correlation
#' among local regression coefficients in geographically weighted regression.
#' Journal of Geographical Systems, 7(2), 161-187.
#'
#' @examples
#' data(LondonHP)
#' diag <- gwr_collin_diag(
#'   PURCHASE ~ FLOORSZ + UNEMPLOY + PROF,
#'   data = LondonHP,
#'   bw = 64,
#'   adaptive = TRUE
#' )
#' print(diag)
#'
#' @importFrom stats model.frame model.extract model.matrix terms cov.wt
#' @export
gwr_collin_diag <- function(
    formula,
    data,
    bw,
    adaptive = FALSE,
    kernel = c("bisquare", "gaussian", "exp", "tricube", "boxcar"),
    longlat = FALSE,
    p = 2.0,
    theta = 0.0,
    dMat = NULL
) {
    kernel <- match.arg(kernel)

    ### Extract coords
    coords <- as.matrix(sf::st_coordinates(sf::st_centroid(data)))
    n <- nrow(coords)

    ### Extract variables
    mf <- model.frame(formula = formula(formula), data = sf::st_drop_geometry(data))
    mt <- attr(mf, "terms")
    y <- model.extract(mf, "response")
    x <- model.matrix(mt, mf)

    has_intercept <- attr(terms(mf), "intercept") == 1
    indep_vars <- colnames(x)
    indep_vars[which(indep_vars == "(Intercept)")] <- "Intercept"
    colnames(x) <- indep_vars

    var_n <- ncol(x)  # number of variables including intercept

    if (var_n <= 1) {
        stop("The number of independent variables must be larger than one")
    }

    # x1: design matrix without intercept (for VIF calculation)
    if (has_intercept) {
        x1 <- x[, -1, drop = FALSE]
    } else {
        x1 <- x
    }

    ### Compute distance matrix if not provided
    if (is.null(dMat)) {
        if (longlat) {
            # Haversine formula for great circle distance
            lon_rad <- coords[, 1] * pi / 180
            lat_rad <- coords[, 2] * pi / 180

            dMat <- matrix(0, n, n)
            for (i in 1:(n - 1)) {
                for (j in (i + 1):n) {
                    dlat <- lat_rad[j] - lat_rad[i]
                    dlon <- lon_rad[j] - lon_rad[i]
                    a <- sin(dlat / 2)^2 + cos(lat_rad[i]) * cos(lat_rad[j]) * sin(dlon / 2)^2
                    dMat[i, j] <- 2 * 6371 * asin(sqrt(a))  # Earth radius 6371 km
                    dMat[j, i] <- dMat[i, j]
                }
            }
        } else {
            # Euclidean or Minkowski distance
            if (p == 2.0 && theta == 0.0) {
                dMat <- as.matrix(stats::dist(coords))
            } else {
                # Minkowski with rotation
                coords_rot <- coords
                if (theta != 0) {
                    cos_t <- cos(theta)
                    sin_t <- sin(theta)
                    coords_rot[, 1] <- coords[, 1] * cos_t - coords[, 2] * sin_t
                    coords_rot[, 2] <- coords[, 1] * sin_t + coords[, 2] * cos_t
                }
                dMat <- as.matrix(stats::dist(coords_rot, method = "minkowski", p = p))
            }
        }
    }

    ### Kernel weight function
    gw_weight <- function(dist, bw, kernel, adaptive) {
        if (adaptive) {
            # Adaptive: bw is number of neighbors
            bw_dist <- sort(dist)[min(bw, length(dist))]
        } else {
            bw_dist <- bw
        }

        w <- switch(kernel,
            "gaussian" = exp(-0.5 * (dist / bw_dist)^2),
            "exp" = exp(-dist / bw_dist),
            "bisquare" = ifelse(dist <= bw_dist, (1 - (dist / bw_dist)^2)^2, 0),
            "tricube" = ifelse(dist <= bw_dist, (1 - (dist / bw_dist)^3)^3, 0),
            "boxcar" = ifelse(dist <= bw_dist, 1, 0)
        )
        return(w)
    }

    ### Initialize output matrices
    # Correlation matrix: (var_n-1)*var_n/2 pairs
    n_corr_pairs <- (var_n - 1) * var_n / 2
    corr_mat <- matrix(0, nrow = n, ncol = n_corr_pairs)

    # VIF matrix: var_n-1 variables (excluding intercept)
    vifs_mat <- matrix(0, nrow = n, ncol = var_n - 1)

    # Condition index: var_n values per location
    vdp_idx <- matrix(0, nrow = n, ncol = var_n)

    # VDP array: n x var_n x var_n
    vdp_pi <- array(0, dim = c(n, var_n, var_n))

    ### Compute diagnostics for each location
    for (i in 1:n) {
        # Get distances to location i
        dist_vi <- dMat[, i]

        # Compute weights
        W_i <- gw_weight(dist_vi, bw, kernel, adaptive)

        # Normalize weights (sum to 1) - as in GWmodel
        sum_w <- sum(W_i)
        Wi <- W_i / sum_w

        ### Local correlations between all variable pairs
        tag <- 0
        for (j in 1:(var_n - 1)) {
            for (k in (j + 1):var_n) {
                tag <- tag + 1
                corr_mat[i, tag] <- cov.wt(cbind(x[, j], x[, k]),
                                           wt = Wi, cor = TRUE)$cor[1, 2]
            }
        }

        ### VIF calculation using weighted correlation matrix (without intercept)
        # This is the standard approach: VIF = diag(solve(R))
        # where R is the correlation matrix
        tryCatch({
            corr_mati <- cov.wt(x1, wt = Wi, cor = TRUE)$cor
            vifs_mat[i, ] <- diag(solve(corr_mati))
        }, error = function(e) {
            vifs_mat[i, ] <- NA
        })

        ### Variance-decomposition proportions (VDP)
        # Following Belsley, Kuh, and Welsch (1980)
        # 1. Weight the design matrix: X * W (element-wise by row)
        xw <- sweep(x, 1, Wi, "*")

        # 2. Column-scale to unit length: divide each column by its L2 norm
        col_norms <- sqrt(colSums(xw^2))
        col_norms[col_norms == 0] <- 1  # avoid division by zero
        xw_scaled <- sweep(xw, 2, col_norms, "/")

        # 3. SVD of scaled weighted design matrix
        svd_x <- svd(xw_scaled)

        # 4. Condition indices: max(d) / d[k]
        vdp_idx[i, ] <- svd_x$d[1] / pmax(svd_x$d, 1e-10)

        # 5. VDP calculation
        # Phi = V * diag(1/d), then square and transpose
        Phi <- svd_x$v %*% diag(1 / pmax(svd_x$d, 1e-10))
        Phi <- t(Phi^2)

        # 6. Normalize by column (prop.table by column)
        col_sums <- colSums(Phi)
        col_sums[col_sums == 0] <- 1
        pi_ij <- sweep(Phi, 2, col_sums, "/")

        vdp_pi[i, , ] <- pi_ij
    }

    ### Create column names for correlations
    corr_nms <- c()
    for (i in 1:(var_n - 1)) {
        for (j in (i + 1):var_n) {
            corr_nms <- c(corr_nms, paste("Corr", paste(indep_vars[i], indep_vars[j], sep = "."), sep = "_"))
        }
    }

    ### Extract final outputs (following GWmodel convention)
    # local_CN: the largest condition index (last column)
    local_CN <- vdp_idx[, var_n]

    # VDP: proportions for the smallest singular value (last column of each pi_ij)
    VDP <- vdp_pi[, var_n, ]

    ### Create result data frame
    # VIF column names (without intercept)
    vif_nms <- paste(indep_vars[-1], "VIF", sep = "_")

    # VDP column names (all variables)
    vdp_nms <- paste(indep_vars, "VDP", sep = "_")

    res_df <- data.frame(
        vifs_mat,
        local_CN = local_CN,
        VDP,
        corr_mat
    )
    colnames(res_df) <- c(vif_nms, "local_CN", vdp_nms, corr_nms)

    ### Create sf object
    res_df$geometry <- sf::st_geometry(data)
    SDF <- sf::st_sf(res_df)

    ### Return results
    result <- list(
        SDF = SDF,
        corr.mat = corr_mat,
        VIF = vifs_mat,
        local_CN = local_CN,
        VDP = VDP,
        indep_vars = indep_vars,
        call = match.call()
    )
    class(result) <- "gwcollindiag"
    result
}

#' Print GWR Collinearity Diagnostics
#' @param x A gwcollindiag object
#' @param ... Additional arguments
#' @method print gwcollindiag
#' @export
print.gwcollindiag <- function(x, ...) {
    cat("GWR Collinearity Diagnostics\n")
    cat("============================\n")
    cat("Formula:", deparse(x$call$formula), "\n")
    cat("Variables:", paste(x$indep_vars, collapse = ", "), "\n\n")

    cat("Local Condition Number (CN) Summary:\n")
    print(summary(x$local_CN))
    cat("\n")

    cat("Local VIF Summary (excluding intercept):\n")
    vif_summary <- apply(x$VIF, 2, summary)
    colnames(vif_summary) <- x$indep_vars[-1]
    print(vif_summary)
    cat("\n")

    n_high_cn <- sum(x$local_CN > 30, na.rm = TRUE)
    n_high_vif <- sum(apply(x$VIF, 1, function(v) any(v > 10, na.rm = TRUE)), na.rm = TRUE)

    cat("Diagnostic Summary:\n")
    cat(sprintf("  Locations with CN > 30: %d (%.1f%%)\n",
                n_high_cn, 100 * n_high_cn / length(x$local_CN)))
    cat(sprintf("  Locations with any VIF > 10: %d (%.1f%%)\n",
                n_high_vif, 100 * n_high_vif / length(x$local_CN)))

    cat("\nLocal VDP Summary (for smallest singular value):\n")
    vdp_summary <- apply(x$VDP, 2, summary)
    colnames(vdp_summary) <- x$indep_vars
    print(vdp_summary)
}

#' Plot GWR Collinearity Diagnostics
#' @param x A gwcollindiag object
#' @param which Which plot: "CN" for condition number, "VIF" for VIF maps, "VDP" for VDP maps
#' @param ... Additional arguments passed to sf::plot
#' @method plot gwcollindiag
#' @export
plot.gwcollindiag <- function(x, which = "CN", ...) {
    if (which == "CN") {
        plot(x$SDF["local_CN"], main = "Local Condition Number", ...)
    } else if (which == "VIF") {
        vif_cols <- grep("_VIF$", names(x$SDF), value = TRUE)
        plot(x$SDF[vif_cols], ...)
    } else if (which == "VDP") {
        vdp_cols <- grep("_VDP$", names(x$SDF), value = TRUE)
        plot(x$SDF[vdp_cols], ...)
    }
}
