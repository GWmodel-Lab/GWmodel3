reg_formula <- function(dep_var, indep_vars) {
    paste(dep_var, paste(indep_vars, collapse = "+"), sep = "~")
}