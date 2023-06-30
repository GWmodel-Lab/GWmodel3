reg_formula <- function(dep_var, indep_vars) {
    paste(dep_var, paste(indep_vars, collapse = "+"), sep = "~")
}

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

enum_list <- function(x, labels, default = labels[[1]]) {
    do.call(switch, c(list(x), labels, default))
}

parallel_types <- list(
    "none" = 1,
    "omp" = 2,
    "cuda" = 4,
    "cluster" = 8
)

kernel_enums <- c(
    "gaussian",
    "exp",
    "bisquare",
    "tricube",
    "boxcar"
)
