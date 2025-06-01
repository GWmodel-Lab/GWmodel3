// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "utils.h"
#include "gwmodel.h"


using namespace std;
using namespace Rcpp;
using namespace arma;
using namespace gwm;

// [[Rcpp::export]]
List gwr_lcr_fit(
    const arma::mat& x,
    const arma::vec& y,
    const arma::mat& coords,
    double bw,
    bool adaptive,
    size_t kernel,
    bool longlat,
    double p,
    double theta,
    double lambda,
    bool lambda_adjust,
    double cn_thresh,
    bool intercept,
    bool hatmatrix,
    size_t parallel_type,
    const IntegerVector& parallel_arg,
    bool optim_bw
) {
    vector<int> vpar_args = as< vector<int> >(IntegerVector(parallel_arg));

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

    GWRLocalCollinearity algorithm;
    algorithm.setCoords(coords);
    algorithm.setDependentVariable(y);
    algorithm.setIndependentVariables(x);
    algorithm.setSpatialWeight(spatial);
    algorithm.setHasIntercept(intercept);
    algorithm.setHasHatMatrix(hatmatrix);

    algorithm.setLambdaAdjust(lambda_adjust);
    algorithm.setLambda(lambda);
    algorithm.setCnThresh(cn_thresh);

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

    if (optim_bw)
    {
        algorithm.setIsAutoselectBandwidth(true);
        algorithm.setBandwidthSelectionCriterion(GWRLocalCollinearity::BandwidthSelectionCriterionType::CV);
    }



    algorithm.fit();
    
    // Return Results
    mat betas = algorithm.betas();
    vec yhat = sum(betas % x, 1);
    List result_list = List::create(
        Named("betas") = betas,
        Named("diagnostic") = mywrap(algorithm.diagnostic()),
        Named("yhat") = yhat,
        Named("localCN") = algorithm.localCN(),
        Named("localLambda") = algorithm.localLambda()
    );
    
    if (optim_bw){
        double bw_value = algorithm.spatialWeight().weight<BandwidthWeight>()->bandwidth();
        result_list["bandwidth"] = wrap(bw_value);
    }

    return result_list;
}
