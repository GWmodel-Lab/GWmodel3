// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "utils.h"
#include "gwmodel.h"

using namespace std;
using namespace Rcpp;
using namespace arma;
using namespace gwm;

// [[Rcpp::export]]
List gwpca_cal(
    const arma::mat& x,
    const arma::mat& coords,
    double bw,
    bool adaptive,
    size_t kernel,
    bool longlat,
    double p,
    double theta,
    size_t keep_components
) {
    // Make Spatial Weight
    BandwidthWeight bandwidth(bw, adaptive, BandwidthWeight::KernelFunctionType((size_t)kernel));
    Distance* distance = nullptr;
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

    // Make Algorithm Object
    GWPCA algorithm;
    algorithm.setCoords(coords);
    algorithm.setVariables(x);
    algorithm.setSpatialWeight(spatial);
    algorithm.setKeepComponents(keep_components);
    try
    {
        algorithm.run();
    }
    catch(const std::exception& e)
    {
        stop(e.what());
    }

    // Return Results
    cube loadings = algorithm.loadings();

    // defined but not calculate in libgwmodel
    // cube scores = algorithm.scores();

    List result_list = List::create(
        Named("local_loadings") = loadings,
        Named("local_PV") = algorithm.localPV()
        // Named("local_scores") = scores
    );
    return result_list;
}