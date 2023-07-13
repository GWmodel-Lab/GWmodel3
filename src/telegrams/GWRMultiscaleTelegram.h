#ifndef GWRMULTISCALETELEGRAM
#define GWRMULTISCALETELEGRAM

#include "RTelegram.h"
#include <map>

class GWRMultiscaleTelegram : public RTelegram
{
public:
    enum class InfoTag {
        Stage,
        BandwidthCriterion,
        InitialBandwidth,
        Backfitting
    };
    static std::map<std::string, InfoTag> TagDict;
    enum class BackfittingInfoTag {
        Iteration,
        VariableBandwidthSelection,
        BackfittingCriterion,
        Stage
    };
    static std::map<std::string, BackfittingInfoTag> BackfittingTagDict;
    static std::map<gwm::GWRMultiscale::BandwidthSelectionCriterionType, std::string> BwCriterionName;

public:
    GWRMultiscaleTelegram(const gwm::GWRMultiscale& algorithm, std::vector<std::string> varNames) : RTelegram(), mAlgorithm(algorithm), mVariableNames(varNames) {}
    ~GWRMultiscaleTelegram() {}
    void parseInfo(std::string message) override;
    void parseBackfittingInfo(std::vector<std::string> messages);
    bool splitInitialBandwidth(const std::string& s, size_t& variable, double& criterion);

private:
    const gwm::GWRMultiscale& mAlgorithm;
    bool mIsInitialBandwidthStage = false;
    bool mIsBackfittingBandwidth = false;
    size_t mCurrentIteration;
    size_t mCurrentVariable;
    std::vector<std::string> mVariableNames;
};

#endif  // GWRMULTISCALETELEGRAM
