// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "utils.h"
#include "gwmodel.h"

using namespace std;
using namespace Rcpp;
using namespace arma;
using namespace gwm;

// [[Rcpp::export]]
List gwr_scalable_fit(
    const arma::mat& x,
    const arma::vec& y,
    const arma::mat& coords,
    double bw,
    bool adaptive,
    size_t kernel,
    bool longlat,
    double p,
    double theta,
    double polynomial,
    bool hatmatrix,
    bool intercept,
    bool optim_bw,
    size_t optim_bw_criterion
)
{
    // Make Spatial Weight
    BandwidthWeight bandwidth(bw, adaptive, BandwidthWeight::KernelFunctionType((size_t)kernel));

    Distance *distance = nullptr;
    if (longlat)
    {
        distance = new CRSDistance(true);
    }
    else
    {
        if (p == 2.0 && theta == 0.0)
        {
            distance = new CRSDistance(false);
        }
        else
        {
            distance = new MinkwoskiDistance(p, theta);
        }
    }

    SpatialWeight spatial(&bandwidth, distance);

    GWRScalable algorithm;
    algorithm.setCoords(coords);
    algorithm.setDependentVariable(y);
    algorithm.setIndependentVariables(x);
    algorithm.setSpatialWeight(spatial);

    algorithm.setPolynomial(arma::uword(polynomial));
    algorithm.setHasHatMatrix(hatmatrix);
    algorithm.setHasIntercept(intercept);
    
    if (optim_bw)
    {
        algorithm.setParameterOptimizeCriterion(GWRScalable::BandwidthSelectionCriterionType(size_t(optim_bw_criterion)));
    }
    try
    {
        algorithm.fit();
    }
    catch(const std::exception& e)
    {
        stop(e.what());
    }

    List results;
    results = List::create(
        Named("diagnostic") = mywrap(algorithm.diagnostic()),
        Named("betas") = algorithm.betas(),
        Named("cv") = algorithm.cv(),
        Named("scale") = algorithm.scale(),
        Named("penalty") = algorithm.penalty());

    return results;
}