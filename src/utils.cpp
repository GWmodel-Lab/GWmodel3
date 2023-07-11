#include "utils.h"
#include <sstream>
#include <vector>
#include <queue>

using namespace std;
using namespace Rcpp;
using namespace arma;
using namespace gwm;

List mywrap(const RegressionDiagnostic& diagnostic)
{
    return List::create(
        Named("RSS") = diagnostic.RSS,
        Named("AIC") = diagnostic.AIC,
        Named("AICc") = diagnostic.AICc,
        Named("ENP") = diagnostic.ENP,
        Named("EDF") = diagnostic.EDF,
        Named("RSquare") = diagnostic.RSquare,
        Named("RSquareAdjust") = diagnostic.RSquareAdjust
    );
}

List mywrap(const VariablesCriterionList& criterion_list)
{
    List model_combinations;
    NumericVector model_criterions;
    for (auto &&item : criterion_list)
    {
        model_combinations.push_back(wrap(item.first));
        model_criterions.push_back(item.second);
    }
    return List::create(
        Named("models") = model_combinations,
        Named("criterions") = model_criterions
    );
}

void RTelegram::print(std::string message, ITelegram::LogLevel level, std::string fun_name, std::string file_name)
{
    switch (level)
    {
    case Logger::LogLevel::LOG_EMERG:
    case Logger::LogLevel::LOG_ALERT:
    case Logger::LogLevel::LOG_CRIT:
    case Logger::LogLevel::LOG_ERR:
        Rcpp::Rcerr << "ERROR: " << message << " [" << fun_name << "]" << " (in " << file_name << ")\n";
        break;
    case Logger::LogLevel::LOG_WARNING:
    case Logger::LogLevel::LOG_NOTICE:
        Rcpp::Rcout << "MSG: " << message << " [" << fun_name << "]" << " (in " << file_name << ")\n";
        break;
    case Logger::LogLevel::LOG_INFO:
        parseInfo(message);
        break;
    case Logger::LogLevel::LOG_DEBUG:
    default:
        break;
    }
}

vector<string> RTelegram::split(const string& s, const char& sep)
{
    istringstream iss(s);
    vector<string> res;
    string buffer;
    while (getline(iss, buffer, sep)) {
        res.push_back(buffer);
    }
    return res;
}

std::map<std::string, GWRBasicTelegram::InfoTag> GWRBasicTelegram::TagDict = {
    make_pair("#stage", InfoTag::Stage),
    make_pair("#bandwidth-criterion", InfoTag::BandwidthCriterion),
    make_pair("#variable-criterion", InfoTag::VariableCriterion)
};

std::map<gwm::GWRBasic::BandwidthSelectionCriterionType, std::string> GWRBasicTelegram::BwCriterionName = {
    make_pair(GWRBasic::BandwidthSelectionCriterionType::AIC, "AIC"),
    make_pair(GWRBasic::BandwidthSelectionCriterionType::CV, "CV")
};

bool GWRBasicTelegram::splitBandwidthCriterion(const string& s, vector<double>& params)
{
    istringstream iss(s);
    string buffer;
    while (getline(iss, buffer, ','))
    {
        if (buffer == "adaptive" || buffer == "fixed" || buffer == "criterion")
            return false;
        else params.push_back(stod(buffer));
    }
    return true;
}

bool GWRBasicTelegram::splitVariableCriterion(const string &s, vector<size_t> &variables, double &criterion)
{
    istringstream iss(s);
    string buffer;
    vector<string> params = split(s, ',');
    if (params[0] == "variables" || params[1] == "criterion")
        return false;
    else
    {
        vector<string> var_ids = split(params[0], '+');
        for (auto &&i : var_ids)
        {
            variables.push_back(stoul(i));
        }
        criterion = stod(params[1]);
        return true;;
    }
}

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
