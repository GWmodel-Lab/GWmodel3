// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "utils.h"
#include "gwmodelpp/GWDA.h"

using namespace std;
using namespace Rcpp;
using namespace arma;
using namespace gwm;

// [[Rcpp::export]]
List gwda_cal(
    const arma::mat &x,
    std::vector<std::string> &y,
    const arma::mat &coords,
    double bw,
    bool adaptive,
    size_t kernel,
    bool longlat,
    double p,
    double theta,
    bool method,
    size_t parallel_type,
    const IntegerVector &parallel_arg)
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

    GWDA algorithm;
    algorithm.setCoords(coords);
    algorithm.setVariables(x);
    algorithm.setGroup(y);
    algorithm.setSpatialWeight(spatial);

    vector<int> vpar_args = as<vector<int>>(IntegerVector(parallel_arg));
    switch (ParallelType(size_t(parallel_type)))
    {
#ifdef _OPENMP
    case ParallelType::OpenMP:
        algorithm.setParallelType(ParallelType::OpenMP);
        algorithm.setOmpThreadNum(vpar_args[0]);
        break;
#endif
    default:
        algorithm.setParallelType(ParallelType::SerialOnly);
        break;
    }

    algorithm.setIsWqda(method);

    try
    {
        algorithm.run();
    }
    catch(const std::exception& e)
    {
        stop(e.what());
    }

    List results;
    results = List::create(
        Named("correctRate") = algorithm.correctRate(),
        Named("group") = algorithm.group(),
        Named("res") = algorithm.res(),
        Named("probs") = algorithm.probs(),
        Named("pmax") = algorithm.pmax(),
        Named("entropy") = algorithm.entropy());

    return results;
}