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
    GWRBasicTelegram(const gwm::GWRBasic& algorithm) : RTelegram(), mAlgorithm(algorithm) {}
    GWRBasicTelegram(const gwm::GWRBasic& algorithm, std::vector<std::string> varNames) : RTelegram(), mAlgorithm(algorithm), mVariableNames(varNames) {}
    ~GWRBasicTelegram() {}
    void parseInfo(std::string message) override;

private:
    const gwm::GWRBasic& mAlgorithm;
    std::vector<std::string> mVariableNames;
};

#endif  // GWRBASICTELEGRAM
