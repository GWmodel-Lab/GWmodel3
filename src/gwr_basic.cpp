#include <Rcpp.h>
#include <armadillo>
#include <vector>
#include <string>
#include "gwmodel.h"
#include "gwmodelpp/GwmVariable.h"
using namespace std;
using namespace arma;
using namespace Rcpp;

// [[Rcpp::plugins("cpp17")]]
// [[Rcpp::export]]
int gwm_gwr_basic(
    const mat& x, const vec& y, const mat& coords,
    const vector<string>& vxfields, const vector<string>& vyfields,
    double bw, bool adaptive, size_t kernel,
    bool longlat, double p, double theta,
    bool hatmatrix, size_t parallel_type, int parallel_arg)
{
    // Make SimpleLayer
    Rcpp::Rcout << "Make SimpleLayer\n";
    mat data = join_rows(y, x);
    vector<string> vfields = vxfields;
    vfields.insert(vfields.begin(), vyfields.begin(), vyfields.begin() + 1);
    CGwmSimpleLayer layer(coords, data, vfields);

    // Make Variables
    Rcpp::Rcout << "Make Variables\n";
    GwmVariable m_dep_var(0, true, vyfields[0]);
    vector<GwmVariable> m_indep_vars;
    for (size_t i = 0; i < vxfields.size(); i++)
    {
        m_indep_vars.push_back(GwmVariable(i + 1, true, vxfields[i]));
    }
    
    // Make Spatial Weight
    Rcpp::Rcout << "Make Spatial Weight\n";
    CGwmBandwidthWeight bandwidth(bw, adaptive, CGwmBandwidthWeight::KernelFunctionType(kernel));
    CGwmDistance* distance = nullptr;
    if (longlat)
    {
        distance = new CGwmCRSDistance(true);
    }
    else
    {
        if (p == 2.0 && theta == 0.0)
        {
            distance = new CGwmCRSDistance(false);
        }
        else
        {
            distance = new CGwmMinkwoskiDistance(p, theta);
        }
    }
    CGwmSpatialWeight spatial(&bandwidth, distance);
    
    // Make Algorithm Object
    Rcpp::Rcout << "Make Algorithm Object\n";
    CGwmGWRBasic algorithm;
    algorithm.setSourceLayer(layer);
    algorithm.setDependentVariable(m_dep_var);
    algorithm.setIndependentVariables(m_indep_vars);
    algorithm.setSpatialWeight(spatial);
    algorithm.setHasHatMatrix(hatmatrix);
    algorithm.run();

    // Calculate Fitted Values
    Rcpp::Rcout << "Calculate Fitted Values\n";
    mat mbetas = algorithm.betas();
    vec mfitted = CGwmGWRBasic::Fitted(x, mbetas);

    return 0;
}