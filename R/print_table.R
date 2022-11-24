#' Print a character matrix as a table.
#'
#' @param x A character matrix.
#' @param col.sep Column seperator. Default to `""`.
#' @param header.sep Header seperator. Default to `"-"`.
#' @param row.begin Character at the beginning of each row.
#' Default to `col.sep`.
#' @param row.end Character at the ending of each row.
#' Default to `col.sep`.
#' @param table.style Name of pre-defined style.
#' Possible values are `"plain"`, `"md"` or `"latex"`. Default to `"plain"`.
#' @param \dots Additional style control arguments.
#'
#' @return No return.
#'
#' @details
#' When `table.style` is specified, `col.sep`, `header.sep`, `row.begin`
#' and `row.end` would not take effects.
#' Because this function will automatically set their values.
#' For each possible value of `table.style`, its corresponding style settings
#' are shown in the following table.
#' \tabular{llll}{
#'                   \tab \strong{\code{plain}} \tab \strong{\code{md}} \tab \strong{\code{latex}} \cr
#' \code{col.sep}    \tab \code{""}             \tab \code{"|"}         \tab \code{"&"}            \cr
#' \code{header.sep} \tab \code{""}             \tab \code{"-"}         \tab \code{""}             \cr
#' \code{row.begin}  \tab \code{""}             \tab \code{"|"}         \tab \code{""}             \cr
#' \code{row.end}    \tab \code{""}             \tab \code{"|"}         \tab \code{"\\\\"}
#' }
#'
#' In this function, characters are right padded by spaces.
#'
#' @seealso [print.gwrm()].
print_table_md <- function(x, col.sep = "", header.sep = "",
                           row.begin = "", row.end = "",
                           table.style = c("plain", "md", "latex"), ...) {
    if (!missing(table.style)) {
        table.style <- match.arg(table.style)
        if (table.style == "md") {
            col.sep <- "|"
            header.sep <- "-"
            row.begin <- "|"
            row.end <- "|"
        } else if (table.style == "latex") {
            col.sep <- "&"
            header.sep <- ""
            row.begin <- ""
            row.end <- "\\\\"
        } else if (table.style == "plain") {
            col.sep <- ""
            header.sep <- ""
            row.begin <- ""
            row.end <- ""
        } else {
           stop("Unknown table.style.")
        }
    }
    if (nchar(header.sep) > 1) {
       stop("Currently only 1 character header.sep is supported.")
    }
    ### Print table
    x.length <- apply(x, c(1, 2), nchar)
    x.length.max <- apply(x.length, 2, max)
    x.fmt <- sprintf("%%%ds", x.length.max)
    for(c in 1:ncol(x)) {
        if(x.length.max[c] > 0)
            cat(ifelse(c == 1, row.begin, col.sep),
                sprintf(x.fmt[c], x[1, c]), "")
    }
    cat(paste0(row.end, "\n"))
    if (nchar(header.sep) > 0) {
        for(c in 1:ncol(x)) {
            if(x.length.max[c] > 0) {
                header.sep.full <- paste(rep("-", x.length.max[c]),
                                         collapse = "")
                cat(ifelse(c == 1, row.begin, col.sep),
                    sprintf(header.sep.full), "")
            }
        }
        cat(paste0(row.end, "\n"))
    }
    for (r in 2:nrow(x)) {
        for (c in 1:ncol(x)) {
            if(x.length.max[c] > 0)
                cat(ifelse(c == 1, row.begin, col.sep),
                    sprintf(x.fmt[c], x[r, c]), "")
        }
        cat(paste0(row.end, "\n"))
    }
}

#' Convert a numeric matrix to character matrix according to a format string.
#'
#' @param m A numeric matrix.
#' @param fmt Format string. Passing to [base::sprintf()].
#'
#' @seealso [base::sprintf()], [print.hgwrm()], [print.summary.hgwrm()].
matrix2char <- function(m, fmt = "%.3f") {
    mc <- NULL
    if ("array" %in% class(m)) {
        mc <- apply(m, seq(length(dim(m))), function(x) { sprintf(fmt, x) })
    } else {
        mc <- sprintf(fmt, m)
    }
    mc
}
