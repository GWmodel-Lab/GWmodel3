#include <Rcpp.h>
#include <armadillo>
#include "utils.h"
#include "gwmodel.h"

using namespace std;
using namespace Rcpp;
using namespace arma;
using namespace gwm;

// [[Rcpp::export]]
List gwss_fit(
    const NumericMatrix& x,
    const NumericMatrix& coords,
    int mode,
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

    GWSS algorithm(mx, mcoords, spatial);
    algorithm.setGWSSMode(GWSS::GWSSMode(mode));
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

    List results;
    switch (mode)
    {
    case 0:
        results = List::create(
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
        break;
    case 1:
        results = List::create(
            Named("corr") = mywrap(algorithm.localCorr()),
            Named("scorr") = mywrap(algorithm.localSCorr())
        );
        break;
    default:
        break;
    }

    return results;
}

// [[Rcpp::export]]
List gwss_average(
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

    GWSS algorithm(mx, mcoords, spatial);
    algorithm.setGWSSMode(GWSS::GWSSMode::Average);
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