#ifndef GWCORRELATIONTELEGRAM
#define GWCORRELATIONTELEGRAM

#include "RTelegram.h"
#include <map>

class GWCorrelationTelegram : public RTelegram
{
public:
    enum class InfoTag {
        Stage,
        BandwidthCriterion,
        InitialBandwidth
    };
    static std::map<std::string, InfoTag> TagDict;
    static std::map<gwm::GWCorrelation::BandwidthSelectionCriterionType, std::string> BwCriterionName;

public:
    GWCorrelationTelegram(gwm::GWCorrelation& algorithm, std::vector<std::string> varNames, int verbose) : 
        RTelegram(verbose), mVariableNames(varNames)
    {
        mBandwidthCriterionType = algorithm.bandwidthSelectionApproach();
        auto& sws = algorithm.spatialWeights();
        mBandwidthType.resize(sws.size());
        std::transform(sws.cbegin(), sws.cend(), mBandwidthType.begin(), [](const gwm::SpatialWeight& sw)
        {
            return sw.weight<gwm::BandwidthWeight>()->adaptive();
        });
    }
    ~GWCorrelationTelegram() {}
    void parseInfo(std::string message) override;
    bool splitInitialBandwidth(const std::string& s, size_t& variable, double& criterion);

private:
    std::vector<gwm::GWCorrelation::BandwidthSelectionCriterionType> mBandwidthCriterionType;
    std::vector<bool> mBandwidthType;
    bool mIsInitialBandwidthStage = false;
    size_t mCurrentIteration;
    size_t mCurrentVariable;
    std::vector<std::string> mVariableNames;
};

#endif  // GWCORRELATIONTELEGRAM
