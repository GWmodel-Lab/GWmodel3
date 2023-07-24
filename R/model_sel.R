#' Generic method for model auto-selection.
#' @param object Geographically Weighted Model
#' @param \dots Other arguments passing to implementation.
#' 
#' @export
step <- function(object, ...) {
    UseMethod("step")
}

#' @rdname step
#' @export
#' @keywords internal
step.default <- function(object, ...) stats::step(object, ...)

#' @describeIn step Create circle view for
#'  model combinations in model selection.
#' @importFrom grDevices rainbow
#' @importFrom graphics legend lines points text
#' @importFrom stats formula terms
#' @export
step_view_circle <- function(object, ...) {
    if (!inherits(object, "modelselcritl")) {
        stop("This function can only be applied on 'modelselcritl' objects.")
    }

    DeVar <- object$dep_var
    InDeVars <- object$indep_vars
    model.list <- object$models

    n <- length(InDeVars)
    cex <- ifelse(n > 10, 10 / n, 1)

    numModels <- length(model.list)
    alpha <- 2*pi/numModels
    cols <- rainbow(n)
    pchs <- rep(c(8,9,10,15,16,17,18,23,24), length.out = n)
    plot(
        x = 0, y = 0,
        xlim = c(-3*n/4, n+8),
        ylim = c(-n+n/5, n/2),
        cex = 2,
        axes = F,
        pch = 22,
        xlab = "",
        ylab = "",
        ...
    )
    
    for (i in seq_len(numModels)) {
        vars <- attr(terms(formula(model.list[[i]])), which = "term.labels")
        nvar <- length(vars)
        p1 <- c(0, 0)
        for (j in seq_len(nvar)) {
            radius <- sqrt(n) * sqrt(j)
            var.idx <- which(InDeVars == vars[j])
            coord <- c(radius*cos((i-1)*alpha), radius*sin((i-1)*alpha))
            lines(
                x = c(p1[1], coord[1]),
                y = c(p1[2], coord[2]),
                col = "grey",
                lwd = cex
            )
            points(
                x = coord[1],
                y = coord[2],
                col = cols[var.idx],
                pch = pchs[var.idx],
                cex = (cex*i/numModels+0.3)
            )
            p1<-coord
        }
        text(
            x = (radius + 0.5) * cos((i - 1) * alpha),
            y = (radius + 0.5) * sin((i - 1) * alpha),
            labels = as.character(i),
            cex = cex*0.6,
            srt = (i-1)*alpha / pi * 180,
            adj = 0
        )
    }
    legend(
        x = "right",
        col = c("black", cols),
        pch = c(22, pchs),
        c(DeVar, InDeVars),
        box.col="white"
    )
}

#' @describeIn step Create scatter plot for
#'  model selection criterion values
#' @importFrom graphics abline
#' @export 
step_view_value <- function(object, ...) {
    if (!inherits(object, "modelselcritl")) {
        stop("This function can only be applied on 'modelselcritl' objects.")
    }

    ruler <- object$criterion_values
    indep_vars <- object$indep_vars
    plot(
        ruler,
        col = "black",
        pch = 20,
        lty = 5,
        type = "b",
        ylab = object$criterion,
        ...
    )
    num_vars <- length(indep_vars)
    for (i in seq_len(num_vars)) {
       abline(v = sum(num_vars:(num_vars - i + 1)), lty = 2)
    }
}

#' @describeIn step Create scatter plot for
#'  differences of model selection criterion values
#' @param ymin The lower boundary of y-axis.
#' @importFrom graphics abline text
#' @export 
step_view_diff <- function(object, ymin = -50, ...) {
    if (!inherits(object, "modelselcritl")) {
        stop("This function can only be applied on 'modelselcritl' objects.")
    }

    ruler <- object$criterion_values
    ruler_diff <- c(0, diff(ruler))
    threshold <- object$threshold
    if (ymin > 0) ymin = -ymin
    if (-threshold < ymin) ymin = -threshold + ymin
    plot(
        ruler_diff,
        col = "black",
        pch = 20,
        ylab = sprintf("Diff(%s)", object$criterion),
        ylim = c(-50, 0),
        ...
    )
    abline(h = -threshold)
    ruler_diff_show <- which(ruler_diff < -threshold & threshold > ymin)
    text(
        x = ruler_diff_show,
        y = ruler_diff[ruler_diff_show],
        labels = sprintf("%d", ruler_diff_show),
        adj = c(0.5, 1.5)
    )
    indep_vars <- object$indep_vars
    num_vars <- length(indep_vars)
    for (i in seq_len(num_vars)) {
       abline(v = sum(num_vars:(num_vars - i + 1)), lty = 2)
    }
}

#' Plot model selection criterions (the circle view).
#' 
#' @param x An object of `modelselcritl` class.
#' @param y Ignored.
#' @param view Name of view.
#' @param ymin The lower bound of y-axis.
#' @param main The main title.
#' @param \dots Additional parameters passing to [plot()].
#' @method plot modelselcritl
#' @name plot
#' 
#' @examples
#' data(LondonHP)
#' m <- gwr_basic(
#'   PURCHASE ~ FLOORSZ + UNEMPLOY + PROF + BATH2 + BEDS2 +
#'       GARAGE1 + TYPEDETCH + TPSEMIDTCH + TYPETRRD + TYPEBNGLW +
#'       BLDPWW1 +BLDPOSTW + BLD60S + BLD70S + BLD80S + CENTHEAT,
#'   LondonHP, bw = "AIC", adaptive = TRUE
#' )
#' m <- step(m, threshold = 10, bw = Inf, optim_bw = "AIC")
#' plot(m$step)
#' plot(m$step, view = "value")
#' plot(m$step, view = "diff")
#' 
#' @importFrom graphics par
#' @export 
plot.modelselcritl <- function(
    x,
    y,
    ...,
    view = c("circle", "value", "diff"),
    ymin,
    main
) {
    if (!inherits(x, "modelselcritl")) {
        stop("This function can only be applied on 'modelselcritl' objects.")
    }
    view <- match.arg(view)
    if (view == "circle") {
        if (missing(main)) {
            op <- par(mai = c(0, 0, 0, 0), omi = c(0, 0, 0, 0))
            step_view_circle(x, ...)
            par(op)
        } else {
            op <- par(mai = c(0, 0, 1, 0), omi = c(0, 0, 0, 0))
            step_view_circle(x, main = main, ...)
            par(op)
        }
    } else if (view == "value") {
        view_args <- list(x, ...)
        if (!missing(main)) {
            view_args$main <- main
        }
        do.call(step_view_value, view_args)
    } else {
        view_args <- list(x, ...)
        if (!missing(ymin)) {
            view_args$ymin <- ymin
        }
        if (!missing(main)) {
            view_args$main <- main
        }
        do.call(step_view_diff, view_args)
    }
}
