#ifndef GWDRTELEGRAM
#define GWDRTELEGRAM

#include "RTelegram.h"
#include <map>

class GWDRTelegram : public RTelegram
{
public:
    enum class InfoTag {
        Stage,
        BandwidthCriterion,
        VariableCriterion,
    };
    static std::map<std::string, InfoTag> TagDict;
    static std::map<gwm::GWDR::BandwidthCriterionType, std::string> BwCriterionName;

public:
    GWDRTelegram(const gwm::GWDR& algorithm, int verbose) : 
        RTelegram(verbose), mAlgorithm(algorithm) {}
    GWDRTelegram(const gwm::GWDR& algorithm, std::vector<std::string> varNames, int verbose) : 
        RTelegram(verbose), mAlgorithm(algorithm), mVariableNames(varNames) {}
    ~GWDRTelegram() {}
    void parseInfo(std::string message) override;
    bool splitBandwidthCriterion(const std::string& s, std::vector<double>& params);

private:
    const gwm::GWDR& mAlgorithm;
    std::vector<std::string> mVariableNames;
};

#endif  // GWDRTELEGRAM
