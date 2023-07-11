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

class RTelegram : public gwm::ITelegram
{

public:
    RTelegram() : ITelegram() {}

    ~RTelegram() {}

    void print(std::string message, ITelegram::LogLevel level, std::string fun_name, std::string file_name) override;

    void progress(std::size_t current, std::size_t total, std::string fun_name, std::string file_name) override
    {
        (void)current;
        (void)total;
        (void)fun_name;
        (void)file_name;
    }

    void progress(double percent, std::string fun_name, std::string file_name) override
    {
        (void)percent;
        (void)fun_name;
        (void)file_name;
    }

    bool stop() override { return false; }

    std::vector<std::string> split(const std::string& s, const char& sep);

    std::string join(std::vector<std::string> ss, std::string delm);

    virtual void parseInfo(std::string message)
    {
        Rcpp::Rcout << "MSG: " << message << "\n";
    }
};

class GWRBasicTelegram : public RTelegram
{
public:
    enum class InfoTag {
        Stage,
        BandwidthCriterion,
        VariableCriterion,
    };

    static std::map<std::string, InfoTag> TagDict;

    static std::map<gwm::GWRBasic::BandwidthSelectionCriterionType, std::string> BwCriterionName;

public:
    GWRBasicTelegram(const gwm::GWRBasic& algorithm, std::vector<std::string> varNames) : RTelegram(), mAlgorithm(algorithm), mVariableNames(varNames) {}

    ~GWRBasicTelegram() {}

    void parseInfo(std::string message) override;

    bool splitBandwidthCriterion(const std::string& s, std::vector<double>& params);

    bool splitVariableCriterion(const std::string& s, std::vector<std::size_t>& variables, double& criterion);

private:
    const gwm::GWRBasic& mAlgorithm;
    std::vector<std::string> mVariableNames;
};
