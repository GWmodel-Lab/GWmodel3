#include "GWRBasicTelegram.h"

using namespace std;
using namespace Rcpp;
using namespace gwm;

std::map<std::string, GWRBasicTelegram::InfoTag> GWRBasicTelegram::TagDict = {
    make_pair("#stage", InfoTag::Stage),
    make_pair("#bandwidth-criterion", InfoTag::BandwidthCriterion),
    make_pair("#variable-criterion", InfoTag::VariableCriterion)
};

std::map<gwm::GWRBasic::BandwidthSelectionCriterionType, std::string> GWRBasicTelegram::BwCriterionName = {
    make_pair(GWRBasic::BandwidthSelectionCriterionType::AIC, "AIC"),
    make_pair(GWRBasic::BandwidthSelectionCriterionType::CV, "CV")
};

void GWRBasicTelegram::parseInfo(string message)
{
    vector<string> msgs = split(message, ' ');
    InfoTag tag = TagDict[msgs[0]];
    switch (tag) {
    case InfoTag::Stage:
        msgs.erase(msgs.begin());
        Rcout << "* " << join(msgs, " ") << "\n";
        break;
    case InfoTag::BandwidthCriterion:
        {
            if (msgs[1].rfind("adaptive", 0) == 0 || msgs[1].rfind("fixed", 0) == 0)
            {
                Rcout << "* Selecting Bandwidth" << "\n";
                if (mVerbose >= 2) Rcout << "** Bandwidth," << BwCriterionName[mBandwidthCriterion] << "\n";
            }
            else
            {
                if (mVerbose >= 2)
                {
                    Rcout << "** " << msgs[1] << "\n";
                }
            }
            break;
        }
    case InfoTag::VariableCriterion:
        {
            vector<size_t> vars;
            double criterions;
            bool isTitle = !splitVariableCriterion(msgs[1], vars, criterions);
            if (isTitle)
            {
                Rcout << "* Selecting Variables" << "\n";
                if (mVerbose >= 2) Rcout << "** Variables,AIC" << "\n";
            }
            else
            {
                if (mVerbose >= 2)
                {
                    vector<string> vars_name(vars.size());
                    transform(vars.begin(), vars.end(), vars_name.begin(), [this](const size_t index)
                    {
                        return this->mVariableNames[index];
                    });
                    Rcout << "** " << join(vars_name, "+") << "," << criterions << "\n";
                }
            }
        }
        break;
    default:
        break;
    }
}
