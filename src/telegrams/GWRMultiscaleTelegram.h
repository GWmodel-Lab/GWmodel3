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
    GWRMultiscaleTelegram(gwm::GWRMultiscale& algorithm, std::vector<std::string> varNames, int verbose) : 
        RTelegram(verbose), mVariableNames(varNames)
    {
        mBandwidthCriterionType = algorithm.bandwidthSelectionApproach();
        mBackfittingCriterionType = algorithm.criterionType();
        auto& sws = algorithm.spatialWeights();
        mBandwidthType.resize(sws.size());
        std::transform(sws.cbegin(), sws.cend(), mBandwidthType.begin(), [](const gwm::SpatialWeight& sw)
        {
            return sw.weight<gwm::BandwidthWeight>()->adaptive();
        });
    }
    ~GWRMultiscaleTelegram() {}
    void parseInfo(std::string message) override;
    void parseBackfittingInfo(std::vector<std::string> messages);
    bool splitInitialBandwidth(const std::string& s, size_t& variable, double& criterion);

private:
    std::vector<gwm::GWRMultiscale::BandwidthSelectionCriterionType> mBandwidthCriterionType;
    std::vector<bool> mBandwidthType;
    gwm::GWRMultiscale::BackFittingCriterionType mBackfittingCriterionType;
    bool mIsInitialBandwidthStage = false;
    bool mIsBackfittingBandwidth = false;
    size_t mCurrentIteration;
    size_t mCurrentVariable;
    std::vector<std::string> mVariableNames;
};

#endif  // GWRMULTISCALETELEGRAM
