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

string RTelegram::join(vector<string> ss, string delm)
{
    if (ss.size() <= 0) return string();
    string res = ss.front();
    for (auto i = ss.begin() + 1; i != ss.end(); i++)
    {
        res += (delm + *i);
    }
    return res;
}

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
            vector<double> bwc;
            bool isTitle = !splitBandwidthCriterion(msgs[1], bwc);
            if (isTitle)
            {
                Rcout << "* Selecting Bandwidth" << "\n";
                Rcout << "** Bandwidth," << BwCriterionName[mAlgorithm.bandwidthSelectionCriterion()] << "\n";
            }
            else
            {
                if (mAlgorithm.spatialWeight().weight<BandwidthWeight>()->adaptive())
                    Rcout << "** " << (int)bwc[0] << "," << bwc[1] << "\n";
                else
                    Rcout << "** " << bwc[0] << "," << bwc[1] << "\n";
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
                Rcout << "** Variables,AIC" << "\n";
            }
            else
            {
                vector<string> vars_name(vars.size());
                transform(vars.begin(), vars.end(), vars_name.begin(), [this](const size_t index)
                {
                    return this->mVariableNames[index];
                });
                Rcout << "** " << join(vars_name, "+") << "," << criterions << "\n";
            }
        }
        break;
    default:
        break;
    }
}
