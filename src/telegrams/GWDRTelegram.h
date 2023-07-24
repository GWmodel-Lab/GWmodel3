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

    GWDRTelegram(const gwm::GWDR& algorithm, std::vector<std::string> varNames, int verbose) : 
        GWDRTelegram(algorithm, verbose)
    {
        mVariableNames = varNames;
    }

    ~GWDRTelegram() {}
    void parseInfo(std::string message) override;
    bool testBandwidthCriterionTitle(const std::string& s);

private:
    // const gwm::GWDR& mAlgorithm;
    std::size_t mDims;
    gwm::GWDR::BandwidthCriterionType mBandwidthCriterionType;
    std::vector<bool> mBandwidthTypes;
    std::vector<std::string> mVariableNames;
};

#endif  // GWDRTELEGRAM
