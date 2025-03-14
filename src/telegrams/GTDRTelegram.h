#ifndef GTDRTELEGRAM
#define GTDRTELEGRAM

#include "RTelegram.h"
#include <map>

class GTDRTelegram : public RTelegram
{
public:
    enum class InfoTag {
        Stage,
        BandwidthCriterion,
        VariableCriterion,
    };
    static std::map<std::string, InfoTag> TagDict;
    static std::map<gwm::GTDR::BandwidthCriterionType, std::string> BwCriterionName;

public:
    GTDRTelegram(const gwm::GTDR& algorithm, int verbose) : 
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

    GTDRTelegram(const gwm::GTDR& algorithm, std::vector<std::string> varNames, int verbose) : 
        GTDRTelegram(algorithm, verbose)
    {
        mVariableNames = varNames;
    }

    ~GTDRTelegram() {}
    void parseInfo(std::string message) override;
    bool testBandwidthCriterionTitle(const std::string& s);

private:
    // const gwm::GTDR& mAlgorithm;
    std::size_t mDims;
    gwm::GTDR::BandwidthCriterionType mBandwidthCriterionType;
    std::vector<bool> mBandwidthTypes;
    std::vector<std::string> mVariableNames;
};

#endif  // GTDRTELEGRAM
