#ifndef RTELEGRAM
#define RTELEGRAM

#include <vector>
#include <string>
#include <Rcpp.h>
#include "gwmodel.h"

class RTelegram : public gwm::ITelegram
{

public:
    RTelegram(int verbose) : mVerbose(verbose) {}

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
    bool splitBandwidthCriterion(const std::string& s, std::vector<double>& params);
    bool splitVariableCriterion(const std::string& s, std::vector<std::size_t>& variables, double& criterion);

protected:
    int mVerbose = 0;
};

#endif  // RTELEGRAM
