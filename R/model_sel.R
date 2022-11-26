#' Generic method for model auto-selection.
#' @param x Geographically Weighted Model
#' @param \dots Other arguments passing to implementation.
#' 
#' @export
model_sel <- function(x, ...) {
    UseMethod("model_sel")
}

#' Create circle view for model selection criterions
#' 
#' @param object An object of `modelselcritl` class
#' 
#' @export 
model_sel_view_circle <- function(object) {
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
    par(mai = rep(0, times = 4))
    plot(
        x = 0, y = 0,
        xlim = c(-3*n/4, n+8),
        ylim = c(-n+n/5, n/2),
        cex = 2,
        axes = F,
        pch = 22,
        xlab = "",
        ylab = ""
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
