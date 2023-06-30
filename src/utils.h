#include <Rcpp.h>
#include <armadillo>
#include "gwmodel.h"
#include "gwmodelpp/Logger.h"

arma::mat myas(const Rcpp::NumericMatrix& rmat);
arma::vec myas(const Rcpp::NumericVector& rvec);
SEXP mywrap(const arma::mat& amat);
SEXP mywrap(const arma::vec& avec);
Rcpp::List mywrap(const gwm::RegressionDiagnostic& diagnostic);
Rcpp::List mywrap(const gwm::VariablesCriterionList& criterion_list);
void r_printer(std::string message, gwm::Logger::LogLevel level, std::string fun_name, std::string file_name);
