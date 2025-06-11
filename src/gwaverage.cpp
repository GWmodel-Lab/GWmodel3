// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "utils.h"
#include "gwmodel.h"

using namespace std;
using namespace Rcpp;
using namespace arma;
using namespace gwm;

// [[Rcpp::export]]
List gwaverage_fit(
    const arma::mat &x,
    const arma::mat &coords,
    double bw,
    bool quantile,
    bool adaptive,
    size_t kernel,
    bool longlat,
    double p,
    double theta,
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

    GWAverage algorithm;
    algorithm.setCoords(coords);
    algorithm.setVariables(x);
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


    algorithm.run();

    List results;
    results = List::create(
        Named("mean") = algorithm.localMean(),
        Named("var") = algorithm.localVar(),
        Named("sdev") = algorithm.localSDev(),
        Named("skew") = algorithm.localSkewness(),
        Named("cv") = algorithm.localCV());
    if (quantile)
    {
        results["median"] = algorithm.localMedian();
        results["iqr"] = algorithm.iqr();
        results["qi"] = algorithm.qi();
    }
    return results;
}