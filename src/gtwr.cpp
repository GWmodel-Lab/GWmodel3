#include <Rcpp.h>
#include <armadillo>
#include "utils.h"
#include "gwmodelpp/GTWR.h"

using namespace std;
using namespace Rcpp;
using namespace arma;
using namespace gwm;

// [[Rcpp::export]]
List gtwr_fit(
    const NumericMatrix& x,
    const NumericVector& y,
    const NumericMatrix& coords,
    const NumericVector& times,
    double bw,
    bool adaptive,
    int kernel,
    double lambda,
    bool longlat,
    double p,
    double theta,
    bool intercept,
    bool hatmatrix,
    size_t parallel_type,
    const IntegerVector& parallel_arg,
    bool optim_bw,
    size_t optim_bw_criterion,
    bool optim_lambda,
    int verbose
) {
    // Convert data types
    arma::mat mx = myas(x);
    arma::vec my = myas(y);
    arma::mat mcoords = myas(coords);
    arma::mat mtimes = myas(times);
    std::vector<int> vpar_args = as< std::vector<int> >(Rcpp::IntegerVector(parallel_arg));

    // Make Spatial Weight
    BandwidthWeight bandwidth(bw, adaptive, BandwidthWeight::KernelFunctionType(kernel));
    CRSDistance *sdist = nullptr;
    if (longlat)
    {
        sdist = new CRSDistance(true);
    }
    else
    {
        if (p == 2.0 && theta == 0.0)
        {
            sdist = new CRSDistance(false);
        }
        else
        {
            sdist = new MinkwoskiDistance(p, theta);
        }
    }
    OneDimDistance tdist;
    CRSSTDistance stdist(sdist, &tdist, lambda);
    SpatialWeight spatial(&bandwidth, &stdist);
    
    // Make Algorithm Object
    GTWR algorithm;
    algorithm.setDependentVariable(my);
    algorithm.setIndependentVariables(mx);
    algorithm.setCoords(mcoords, mtimes);
    algorithm.setSpatialWeight(spatial);
    algorithm.setHasHatMatrix(hatmatrix);
    algorithm.setHasIntercept(intercept);

    if (optim_bw) 
    {
        algorithm.setIsAutoselectBandwidth(true);
        algorithm.setBandwidthSelectionCriterion(GTWR::BandwidthSelectionCriterionType(size_t(optim_bw_criterion)));
    }

    if (optim_lambda)
    {
        algorithm.setIsAutoselectLambda(true);
    }

    switch (ParallelType(size_t(parallel_type)))
    {
    case ParallelType::SerialOnly:
        algorithm.setParallelType(ParallelType::SerialOnly);
        break;
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
    try
    {
        algorithm.fit();
    }
    catch(const std::exception& e)
    {
        stop(e.what());
    }
        
    // Return Results
    mat betas = algorithm.betas();
    vec fitted = sum(mx % betas, 1);
    List result_list = List::create(
        Named("betas") = mywrap(betas),
        Named("betasSE") = mywrap(algorithm.betasSE()),
        Named("sTrace") = mywrap(algorithm.sHat()),
        Named("sHat") = mywrap(algorithm.s()),
        Named("diagnostic") = mywrap(algorithm.diagnostic()),
        Named("fitted") = mywrap(fitted)
    );
    const SpatialWeight& spatialWeights = algorithm.spatialWeight();
    if (optim_bw)
    {
        result_list["bw_value"] = wrap(spatialWeights.weight<BandwidthWeight>()->bandwidth());
    }
    // if (optim_lambda)
    // {
    //     result_list["lambda_value"] = wrap(spatialWeights.distance<CRSSTDistance>()->mLambda);
    // }

    return result_list;
}
