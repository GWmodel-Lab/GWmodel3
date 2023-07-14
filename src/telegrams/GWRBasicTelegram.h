#ifndef GWRBASICTELEGRAM
#define GWRBASICTELEGRAM

#include "RTelegram.h"
#include <map>

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
    GWRBasicTelegram(gwm::GWRBasic& algorithm, int verbose) : 
        RTelegram(verbose)
    {
        mBandwidthCriterion = algorithm.bandwidthSelectionCriterion();
    }

    GWRBasicTelegram(gwm::GWRBasic& algorithm, std::vector<std::string> varNames, int verbose) : 
        GWRBasicTelegram(algorithm, verbose)
    {
        mVariableNames = varNames;
    }

    ~GWRBasicTelegram() {}
    void parseInfo(std::string message) override;

private:
    gwm::GWRBasic::BandwidthSelectionCriterionType mBandwidthCriterion = gwm::GWRBasic::BandwidthSelectionCriterionType::AIC;
    std::vector<std::string> mVariableNames;
};

#endif  // GWRBASICTELEGRAM
