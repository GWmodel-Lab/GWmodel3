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


// [[Rcpp::export]]
List gw_correlation(
    const NumericMatrix& x,
    const NumericMatrix& y,
    const NumericMatrix& coords,
    const NumericVector& bw,
    const LogicalVector& adaptive,
    const IntegerVector& kernel,
    const LogicalVector& longlat,
    const NumericVector& p,
    const NumericVector& theta,
    const IntegerVector& initial_type,
    const IntegerVector& optim_bw_criterion,
    size_t parallel_type,
    const IntegerVector& parallel_arg
) {
    // Convert data types
    mat mx = myas(x);
    mat my = myas(y);
    mat mcoords = myas(coords);

    // Make Spatial Weight
    size_t nVar = (size_t)mx.n_cols * (size_t)my.n_cols;
    auto vbw = as< vector<double> >(NumericVector(bw));
    auto vadaptive = as< vector<bool> >(LogicalVector(adaptive));
    auto vkernel = as< vector<int> >(IntegerVector(kernel));
    auto vlonglat = as< vector<bool> >(LogicalVector(longlat));
    auto vp = as< vector<double> >(NumericVector(p));
    auto vtheta = as< vector<double> >(NumericVector(theta));
    auto vinitial_type = as< vector<int> >(IntegerVector(initial_type));
    auto voptim_bw_criterion = as< vector<int> >(IntegerVector(optim_bw_criterion));
    vector<SpatialWeight> spatials;
    for (size_t i = 0; i < vbw.size(); i++)
    {
        BandwidthWeight bandwidth(vbw[i], vadaptive[i], BandwidthWeight::KernelFunctionType(vkernel[i]));
        Distance* distance = nullptr;
        if (vlonglat[i]) {
            distance = new CRSDistance(true);
        }
        else
        {
            if (vp[i] == 2.0 && vtheta[i] == 0.0) distance = new CRSDistance(false);
            else distance = new MinkwoskiDistance(vp[i], vtheta[i]);
        }
        spatials.push_back(SpatialWeight(&bandwidth, distance));
    }
    vector<GWCorrelation::BandwidthInitilizeType> bandwidthInitialize(vinitial_type.size());
    transform(vinitial_type.begin(), vinitial_type.end(), bandwidthInitialize.begin(), [](int x) {
        return GWCorrelation::BandwidthInitilizeType(x);
    });
    vector<GWCorrelation::BandwidthSelectionCriterionType> bandwidthSelectionApproach(voptim_bw_criterion.size());
    transform(voptim_bw_criterion.begin(), voptim_bw_criterion.end(), bandwidthSelectionApproach.begin(), [](int x) {
        return GWCorrelation::BandwidthSelectionCriterionType(x);
    });

    GWCorrelation algorithm;
    algorithm.setVariables1(mx);
    algorithm.setVariables2(my);
    algorithm.setCoords(mcoords);
    algorithm.setSpatialWeights(spatials);
    algorithm.setBandwidthInitilize(bandwidthInitialize);
    algorithm.setBandwidthSelectionApproach(bandwidthSelectionApproach);

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

    vector<double> bw_value;
    const vector<SpatialWeight>& spatialWeights = algorithm.spatialWeights();
    for (size_t i = 0; i < nVar; i++)
    {
        bw_value.push_back(spatialWeights[i].weight<BandwidthWeight>()->bandwidth());
    }

    List results = List::create(
        Named("Cov") = mywrap(algorithm.localCov()),
        Named("Corr") = mywrap(algorithm.localCorr()),
        Named("SCorr") = mywrap(algorithm.localSCorr()),
        Named("bw_value") = wrap(bw_value),
    );

    return results;
}