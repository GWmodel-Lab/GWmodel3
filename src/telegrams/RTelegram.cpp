#include "RTelegram.h"

using namespace std;
using namespace Rcpp;
using namespace gwm;

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

bool RTelegram::splitBandwidthCriterion(const string& s, vector<double>& params)
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

bool RTelegram::splitVariableCriterion(const string &s, vector<size_t> &variables, double &criterion)
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
