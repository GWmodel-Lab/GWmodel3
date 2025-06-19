#include "GWCorrelationTelegram.h"

using namespace std;
using namespace Rcpp;
using namespace gwm;

std::map<std::string, GWCorrelationTelegram::InfoTag> GWCorrelationTelegram::TagDict = {
    make_pair("#stage", InfoTag::Stage),
    make_pair("#bandwidth-criterion", InfoTag::BandwidthCriterion),
    make_pair("#initial-bandwidth", InfoTag::InitialBandwidth),
};

std::map<gwm::GWCorrelation::BandwidthSelectionCriterionType, std::string> GWCorrelationTelegram::BwCriterionName = {
    make_pair(GWCorrelation::BandwidthSelectionCriterionType::AIC, "AIC"),
    make_pair(GWCorrelation::BandwidthSelectionCriterionType::CV, "CV")
};

void GWCorrelationTelegram::parseInfo(std::string message)
{
    vector<string> msgs = split(message, ' ');
    InfoTag tag = TagDict[msgs[0]];
    switch (tag) {
    case InfoTag::Stage:
    {
        msgs.erase(msgs.begin());
        Rcout << "* " << join(msgs, " ") << "\n";
        break;
    }
    case InfoTag::InitialBandwidth:
    {
        size_t variable;
        double bw;
        mIsInitialBandwidthStage = true;
        if (!splitInitialBandwidth(msgs[1], variable, bw))
        {
            mCurrentVariable = variable;
            if (mVerbose >= 2) Rcout << "** Now selecting bandwidth for variable " << mVariableNames[mCurrentVariable] << "\n";
        }
        else
        {
            if (mVerbose >= 2) Rcout << "** Bandwidth selected for variable " << mVariableNames[mCurrentVariable] << " is " << bw << "\n";
            // if (mVerbose >= 3) {
                // Rcout << "*** Bandwidth," << BwCriterionName[mBandwidthCriterionType[mCurrentVariable]] << "\n";
                // Rcout << "*** " << msgs[1] << "\n";
            // }
        }
        break;
    }
    case InfoTag::BandwidthCriterion:
    {
        // mIsInitialBandwidthStage = false;
        if (msgs[1].rfind("adaptive", 0) == 0 || msgs[1].rfind("fixed", 0) == 0)
        {
            if (!mIsInitialBandwidthStage)
            {
                Rcout << "* Selecting Bandwidth" << "\n";
                // Rcout << "-----1\n";
            }
            else
            {
                if (mVerbose >= 3) Rcout << "*** " << msgs[1] << "\n";
                // Rcout << "-----2\n";
            }
        }
        else
        {
            if (mVerbose >= 3) Rcout << "** " << msgs[1] << "\n";
            // Rcout << "-----3\n";
        }
        break;
    }
    default:
        break;
    }
}

bool GWCorrelationTelegram::splitInitialBandwidth(const std::string &s, size_t &variable, double &criterion)
{
    auto parts = split(s, ',');
    variable = stoul(parts[0]);
    if (parts.size() > 1)
    {
        criterion = stod(parts[1]);
        return true;
    }
    else return false;
}
