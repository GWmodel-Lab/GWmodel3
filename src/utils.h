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

inline arma::mat myas(const Rcpp::NumericMatrix& rmat)
{
    return arma::mat(rmat.begin(), rmat.nrow(), rmat.ncol());
}

inline arma::vec myas(const Rcpp::NumericVector& rvec)
{
    return arma::vec(rvec.begin(), rvec.size());
}

inline SEXP mywrap(const arma::mat& amat)
{
    Rcpp::RObject x = Rcpp::wrap(amat.begin(), amat.end());
    x.attr("dim") = Rcpp::Dimension(amat.n_rows, amat.n_cols);
    return x;
}

inline SEXP mywrap(const arma::vec& avec)
{
    Rcpp::RObject x = Rcpp::wrap(avec.begin(), avec.end());
    return x;
}
