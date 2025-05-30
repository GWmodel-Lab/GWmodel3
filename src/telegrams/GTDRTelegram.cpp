#include "GTDRTelegram.h"

using namespace std;
using namespace Rcpp;
using namespace gwm;

map<string, GTDRTelegram::InfoTag> GTDRTelegram::TagDict = {
    make_pair("#stage", InfoTag::Stage),
    make_pair("#bandwidth-criterion", InfoTag::BandwidthCriterion),
    make_pair("#variable-criterion", InfoTag::VariableCriterion)
};

map<GTDR::BandwidthCriterionType, string> GTDRTelegram::BwCriterionName = {
    make_pair(GTDR::BandwidthCriterionType::AIC, "AIC"),
    make_pair(GTDR::BandwidthCriterionType::CV, "CV")
};

void GTDRTelegram::parseInfo(string message)
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
            if (testBandwidthCriterionTitle(msgs[1]))
            {
                Rcout << "* Selecting Bandwidth" << "\n";
                if (mVerbose >= 2)
                {
                    vector<string> titles;
                    for (size_t i = 1; i <= mDims; i++)
                    {
                        titles.push_back(string("Bandwidth") + to_string(i));
                    }
                    titles.push_back(BwCriterionName[mBandwidthCriterionType]);
                    Rcout << "** " << join(titles, ",") << "\n";
                }
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

bool GTDRTelegram::testBandwidthCriterionTitle(const std::string &s)
{
    istringstream iss(s);
    string buffer;
    while (getline(iss, buffer, ','))
    {
        if (buffer.find(':', 0) != string::npos)
            return true;
    }
    return false;
}
