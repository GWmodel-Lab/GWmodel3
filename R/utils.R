reg_formula <- function(dep_var, indep_vars) {
    paste(dep_var, paste(indep_vars, collapse = "+"), sep = "~")
}

.enum <- function(x, labels, values = (seq_along(labels) - 1)) {
    values[match(x, labels, nomatch = 1)]
}
