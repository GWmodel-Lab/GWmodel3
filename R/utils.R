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
#' @rdname print
#' @export 
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
#' @rdname print
#' @noRd 
matrix2char <- function(m, fmt = "%.3f") {
    mc <- NULL
    if ("array" %in% class(m)) {
        mc <- apply(m, seq(length(dim(m))), function(x) { sprintf(fmt, x) })
    } else {
        mc <- sprintf(fmt, m)
    }
    mc
}

#' Create formula by dependent and independent variables
#' 
#' @param dev_var The name of the dependent variable.
#' @param indep_vars Names of independent variables.
#' 
#' @return A formula according to names of variables
#' 
#' @examples
#' reg_formula("y", c("x1", "x2", "x3"))
#' 
#' @noRd 
#' 
reg_formula <- function(dep_var, indep_vars) {
    paste(dep_var, paste(indep_vars, collapse = "+"), sep = "~")
}

#' Convert to enum names to values
#' 
#' @param x The enum name to be converted.
#' @param labels All possible enum names.
#' @param values Enum values corresponding to names.
#' 
#' @return The corresponding enum value.
#' If the name cannot be matched to any one of the labels,
#' return the first value.
#' 
#' @examples 
#' enum("a", c("a", "b", "c"))
#' 
#' enum("a", c("a", "b", "c"), c(1, 2, 4))
#' 
#' enum("e", c("a", "b", "c"), c(1, 2, 4))
#' 
#' @noRd 
#' 
enum <- function(x, labels, values = (seq_along(labels) - 1)) {
    if (missing(labels)) {
        formal.labels <- formals(sys.function(sysP <- sys.parent()))
        labels <- eval(
            formal.labels[[as.character(substitute(x))]],
            envir = sys.frame(sysP)
        )
    }
    values[match(x, labels, nomatch = 1)]
}

#' Convert to enum names to values
#' 
#' @param x The enum name to be converted.
#' @param labels All possible enum names.
#' @param default Default value if \code{x} is missing
#' in names of \code{labels}.
#' 
#' @return The corresponding enum value.
#' If the name cannot be matched to any one of the labels,
#' return the default value.
#' 
#' @examples 
#' enum("a", list("a" = 1, "b" = 2))
#' 
#' enum("e", list("a" = 1, "b" = 2), 1)
#' 
#' @noRd
#' 
enum_list <- function(x, labels, default = labels[[1]]) {
    do.call(switch, c(list(x), labels, default))
}

parallel_types <- list(
    "none" = 1,
    "omp" = 2,
    "cuda" = 4,
    "cluster" = 8
)

kernel_enums <- c("gaussian", "exp", "bisquare", "tricube", "boxcar")
