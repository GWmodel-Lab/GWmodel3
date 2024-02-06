#ifndef SDRTELEGRAM
#define SDRTELEGRAM

#include "RTelegram.h"
#include <map>

class SDRTelegram : public RTelegram
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
    SDRTelegram(const gwm::GWDR& algorithm, int verbose) : 
        RTelegram(verbose)
    {
        const auto& sws = algorithm.spatialWeights();
        mDims = sws.size();
        mBandwidthTypes.resize(mDims);
        std::transform(sws.cbegin(), sws.cend(), mBandwidthTypes.begin(), [](const gwm::SpatialWeight& sw)
        {
            return sw.weight<gwm::BandwidthWeight>()->adaptive();
        });
        mBandwidthCriterionType = algorithm.bandwidthCriterionType();
    }

    SDRTelegram(const gwm::GWDR& algorithm, std::vector<std::string> varNames, int verbose) : 
        SDRTelegram(algorithm, verbose)
    {
        mVariableNames = varNames;
    }

    ~SDRTelegram() {}
    void parseInfo(std::string message) override;
    bool testBandwidthCriterionTitle(const std::string& s);

private:
    // const gwm::GWDR& mAlgorithm;
    std::size_t mDims;
    gwm::GWDR::BandwidthCriterionType mBandwidthCriterionType;
    std::vector<bool> mBandwidthTypes;
    std::vector<std::string> mVariableNames;
};

#endif  // SDRTELEGRAM
