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
