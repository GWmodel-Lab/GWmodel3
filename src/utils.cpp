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
