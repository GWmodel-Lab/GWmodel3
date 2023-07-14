#include "GWRMultiscaleTelegram.h"

using namespace std;
using namespace Rcpp;
using namespace gwm;

std::map<std::string, GWRMultiscaleTelegram::InfoTag> GWRMultiscaleTelegram::TagDict = {
    make_pair("#stage", InfoTag::Stage),
    make_pair("#bandwidth-criterion", InfoTag::BandwidthCriterion),
    make_pair("#initial-bandwidth", InfoTag::InitialBandwidth),
    make_pair("#back-fitting", InfoTag::Backfitting)
};

std::map<std::string, GWRMultiscaleTelegram::BackfittingInfoTag> GWRMultiscaleTelegram::BackfittingTagDict = {
    make_pair("#iteration", BackfittingInfoTag::Iteration),
    make_pair("#variable-bandwidth-selection", BackfittingInfoTag::VariableBandwidthSelection),
    make_pair("#backfitting-criterion", BackfittingInfoTag::BackfittingCriterion)
};

std::map<gwm::GWRMultiscale::BandwidthSelectionCriterionType, std::string> GWRMultiscaleTelegram::BwCriterionName = {
    make_pair(GWRMultiscale::BandwidthSelectionCriterionType::AIC, "AIC"),
    make_pair(GWRMultiscale::BandwidthSelectionCriterionType::CV, "CV")
};

void GWRMultiscaleTelegram::parseInfo(std::string message)
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
        }
        break;
    }
    case InfoTag::BandwidthCriterion:
    {
        if (msgs[1].rfind("adaptive", 0) == 0 || msgs[1].rfind("fixed", 0) == 0)
        {
            if (!mIsInitialBandwidthStage)
            {
                Rcout << "* Selecting Bandwidth" << "\n";
            }
            if (mIsBackfittingBandwidth)
            {
                if (mVerbose >= 4) Rcout << "**** " << "Bandwidth," << BwCriterionName[mBandwidthCriterionType[mCurrentVariable]] << "\n";
            }
            else
            {
                if (mVerbose >= 3) Rcout << "*** " << "Bandwidth," << BwCriterionName[mBandwidthCriterionType[mCurrentVariable]] << "\n";
            }
        }
        else
        {
            if (mIsInitialBandwidthStage)
            {
                if (mVerbose >= 3) Rcout << "*** " << msgs[1] << "\n";
            }
            else if (mIsBackfittingBandwidth)
            {
                if (mVerbose >= 4) Rcout << "**** " << msgs[1] << "\n";
            }
            else
            {
                if (mVerbose >= 2) Rcout << "** " << msgs[1] << "\n";
            }
        }
        break;
    }
    case InfoTag::Backfitting:
    {
        mIsInitialBandwidthStage = false;
        msgs.erase(msgs.begin());
        parseBackfittingInfo(msgs);
        break;
    }
    default:
        break;
    }
}

void GWRMultiscaleTelegram::parseBackfittingInfo(std::vector<std::string> messages)
{
    BackfittingInfoTag tag = (messages[0].rfind("#", 0) == 0) ? BackfittingTagDict[messages[0]] : BackfittingInfoTag::Stage;
    switch (tag)
    {
    case BackfittingInfoTag::Stage:
    {
        Rcout << "* " << join(messages, " ") << "\n";
        break;
    }
    case BackfittingInfoTag::Iteration:
    {
        mCurrentIteration = stoul(messages[1]);
        if (mVerbose >= 2) Rcout << "** Iteration " << mCurrentIteration << "\n";
        break;
    }
    case BackfittingInfoTag::VariableBandwidthSelection:
    {
        vector<string> args = split(messages[1], ',');
        mCurrentVariable = stoul(args[0]);
        mIsBackfittingBandwidth = true;
        string variable = mVariableNames[mCurrentVariable];
        int bwSizePrecision = mBandwidthType[mCurrentVariable] ? 0 : 6;
        switch (args.size())
        {
        case 1:
            if (mVerbose >= 3) Rcout << "*** Now select an optimum bandwidth for the variable " << variable << "\n";
            if (mVerbose >= 4) Rcout << "**** Bandwidth," << BwCriterionName[mBandwidthCriterionType[mCurrentVariable]] << "\n";
            break;
        case 5:
        {
            double bwi0s = stod(args[1]), bwi1s = stod(args[2]), dbw = stod(args[3]);
            bool converged = (args[4] == "true");
            if (mVerbose >= 3)
            {
                Rcout << "*** The newly selected bandwidth for variable " << variable 
                    << " is " << setprecision(bwSizePrecision) << bwi1s 
                    << " (last is " << bwi0s 
                    << ", difference is " << setprecision(6) << dbw << ")\n";
                if (!converged)
                    Rcout << "*** The bandwidth for variable " << variable << " will be continually selected in the text iteration\n";
                else
                    Rcout << "*** The bandwidth for variable " << variable << " seems to be converged and will be kept the same in the following iterations\n";
            }
            break;
        }
        case 6:
        {
            size_t times = stoul(args[4]), retry = stoul(args[5]);
            if (mVerbose >= 3) 
                Rcout << "*** The bandwidth for variable " << variable << " seems to be converted for " << times << " times"
                      << "It will be continually optimized in the next " << retry << " times\n";
            break;
        }
        default:
            break;
        }
        break;
    }
    case BackfittingInfoTag::BackfittingCriterion:
    {
        double criterionValue = stod(messages[1]);
        string criterionName = (mBackfittingCriterionType == GWRMultiscale::BackFittingCriterionType::CVR) ? "change value of RSS (CVR)" : "differential change value of RSS (dCVR)";
        if (mVerbose >= 3) Rcout << "*** The " << criterionName << " is " << criterionValue << "\n";
        if (mVerbose >= 3) Rcout << "** End of iteration " << mCurrentIteration << "\n";
        break;
    }
    default:
        break;
    }
}

bool GWRMultiscaleTelegram::splitInitialBandwidth(const std::string &s, size_t &variable, double &criterion)
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
