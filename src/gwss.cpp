#include <Rcpp.h>
#include <armadillo>
#include "utils.h"
#include "gwmodel.h"

using namespace std;
using namespace Rcpp;
using namespace arma;
using namespace gwm;

// [[Rcpp::export]]
List gw_average(
    const NumericMatrix& x,
    const NumericMatrix& coords,
    bool quantile,
    double bw,
    bool adaptive,
    size_t kernel,
    bool longlat,
    double p,
    double theta,
    size_t parallel_type,
    const IntegerVector& parallel_arg
) {
    // Convert data types
    mat mx = myas(x);
    mat mcoords = myas(coords);

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

    GWAverage algorithm(mx, mcoords, spatial);
    algorithm.setQuantile(quantile);

    vector<int> vpar_args = as< vector<int> >(IntegerVector(parallel_arg));
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

    List results = List::create(
        Named("mean") = mywrap(algorithm.localMean()),
        Named("var") = mywrap(algorithm.localVar()),
        Named("sdev") = mywrap(algorithm.localSDev()),
        Named("skew") = mywrap(algorithm.localSkewness()),
        Named("cv") = mywrap(algorithm.localCV())
    );
    if (quantile) {
        results["median"] = mywrap(algorithm.localMedian());
        results["iqr"] = mywrap(algorithm.iqr());
        results["qi"] = mywrap(algorithm.qi());
    }

    return results;
}