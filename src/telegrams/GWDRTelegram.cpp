#include "GWDRTelegram.h"

using namespace std;
using namespace Rcpp;
using namespace gwm;

map<string, GWDRTelegram::InfoTag> GWDRTelegram::TagDict = {
    make_pair("#stage", InfoTag::Stage),
    make_pair("#bandwidth-criterion", InfoTag::BandwidthCriterion),
    make_pair("#variable-criterion", InfoTag::VariableCriterion)
};

map<GWDR::BandwidthCriterionType, string> GWDRTelegram::BwCriterionName = {
    make_pair(GWDR::BandwidthCriterionType::AIC, "AIC"),
    make_pair(GWDR::BandwidthCriterionType::CV, "CV")
};

void GWDRTelegram::parseInfo(string message)
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
            size_t dims = mAlgorithm.spatialWeights().size();
            bool isTitle = !splitBandwidthCriterion(msgs[1], bwc);
            if (isTitle)
            {
                Rcout << "* Selecting Bandwidth" << "\n";
                if (mVerbose >= 2)
                {
                    vector<string> titles;
                    for (size_t i = 1; i <= dims; i++)
                    {
                        titles.push_back(string("Bandwidth") + to_string(i));
                    }
                    titles.push_back(BwCriterionName[mAlgorithm.bandwidthCriterionType()]);
                    Rcout << "** " << join(titles, ",") << "\n";
                }
            }
            else
            {
                if (mVerbose >= 2)
                {
                    vector<string> values;
                    for (size_t i = 0; i < dims; i++)
                    {
                        string v = mAlgorithm.spatialWeights()[i].weight<BandwidthWeight>()->adaptive() ? to_string((int)bwc[i]) : to_string(bwc[i]);
                        values.push_back(v);
                    }
                    values.push_back(to_string(bwc[dims]));
                    Rcout << "** " << join(values, ",") << "\n";
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
                vector<string> vars_name(vars.size());
                transform(vars.begin(), vars.end(), vars_name.begin(), [this](const size_t index)
                {
                    return this->mVariableNames[index];
                });
                if (mVerbose >= 2) Rcout << "** " << join(vars_name, "+") << "," << criterions << "\n";
            }
        }
        break;
    default:
        break;
    }
}

bool GWDRTelegram::splitBandwidthCriterion(const std::string &s, std::vector<double> &params)
{
    istringstream iss(s);
    string buffer;
    while (getline(iss, buffer, ','))
    {
        if (buffer.find(':', 0) != string::npos)
            return false;
        else params.push_back(stod(buffer));
    }
    return true;
}
