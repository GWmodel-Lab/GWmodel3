// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "gwmodel.h"
#include "gwmodelpp/Logger.h"

Rcpp::List mywrap(const gwm::RegressionDiagnostic& diagnostic);
Rcpp::List mywrap(const gwm::VariablesCriterionList& criterion_list);
