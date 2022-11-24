#include <Rcpp.h>
#include <armadillo>
#include <vector>
#include <string>
#include "gwmodel.h"
using namespace std;
using namespace arma;
using namespace Rcpp;

// [[Rcpp::plugins("cpp17")]]
// [[Rcpp::export]]
int gwm_gwr_basic(
    const mat& x, const vec& y, const mat& coords,
    double bw, bool adaptive, size_t kernel,
    bool longlat, double p, double theta,
    bool hatmatrix, bool intercept, size_t parallel_type, int parallel_arg,
    mat& betas, mat& betasSE, vec& sTrace, mat& sHat, vec& fitted, GwmRegressionDiagnostic& diagnostic)
{    
    // Make Spatial Weight
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
    CGwmGWRBasic algorithm(x, y, coords, spatial, hatmatrix, intercept);
    algorithm.fit();

    // Reture results
    betas = algorithm.betas();
    betasSE = algorithm.betasSE();
    sTrace = algorithm.sHat();
    sHat = algorithm.s();
    fitted = CGwmGWRBasic::Fitted(x, betas);
    diagnostic = algorithm.diagnostic();

    return 0;
}