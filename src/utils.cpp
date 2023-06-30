#include "utils.h"

using namespace Rcpp;
using namespace arma;
using namespace gwm;

inline arma::mat myas(const NumericMatrix& rmat)
{
    return arma::mat(rmat.begin(), rmat.nrow(), rmat.ncol());
}

inline arma::vec myas(const NumericVector& rvec)
{
    return arma::vec(rvec.begin(), rvec.size());
}

inline SEXP mywrap(const arma::mat& amat)
{
    RObject x = wrap(amat.begin(), amat.end());
    x.attr("dim") = Dimension(amat.n_rows, amat.n_cols);
    return x;
}

inline SEXP mywrap(const arma::vec& avec)
{
    RObject x = wrap(avec.begin(), avec.end());
    return x;
}

List mywrap(const RegressionDiagnostic& diagnostic)
{
    return List::create(
        Named("RSS") = diagnostic.RSS,
        Named("AIC") = diagnostic.AIC,
        Named("AICc") = diagnostic.AICc,
        Named("ENP") = diagnostic.ENP,
        Named("EDF") = diagnostic.EDF,
        Named("RSquare") = diagnostic.RSquare,
        Named("RSquareAdjust") = diagnostic.RSquareAdjust
    );
}

List mywrap(const VariablesCriterionList& criterion_list)
{
    List model_combinations;
    NumericVector model_criterions;
    for (auto &&item : criterion_list)
    {
        model_combinations.push_back(wrap(item.first));
        model_criterions.push_back(item.second);
    }
    return List::create(
        Named("models") = model_combinations,
        Named("criterions") = model_criterions
    );
}

void r_printer(std::string message, Logger::LogLevel level, std::string fun_name, std::string file_name)
{
    switch (level)
    {
    case Logger::LogLevel::LOG_EMERG:
    case Logger::LogLevel::LOG_ALERT:
    case Logger::LogLevel::LOG_CRIT:
    case Logger::LogLevel::LOG_ERR:
        Rcpp::Rcerr << "ERROR: " << message << " [" << fun_name << "]" << " (in " << file_name << ")\n";
        break;
    case Logger::LogLevel::LOG_WARNING:
    case Logger::LogLevel::LOG_NOTICE:
    case Logger::LogLevel::LOG_INFO:
    case Logger::LogLevel::LOG_DEBUG:
    default:
        Rcpp::Rcout << "MSG: " << message << " [" << fun_name << "]" << " (in " << file_name << ")\n";
        break;
    }
}
