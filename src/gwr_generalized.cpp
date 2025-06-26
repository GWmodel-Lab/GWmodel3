// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "utils.h"
#include "gwmodel.h"

using namespace std;
using namespace Rcpp;
using namespace arma;
using namespace gwm;

// [[Rcpp::export]]
List gwr_generalized_fit(
    const arma::mat& x,
    const arma::vec& y,
    const arma::mat& coords,
    size_t family,
    double bw,
    bool adaptive,
    size_t kernel,
    bool longlat,
    double p,
    double theta,
    bool hatmatrix,
    bool intercept,
    bool optim_bw,
    size_t optim_bw_criterion,
    size_t parallel_type,
    const IntegerVector& parallel_arg
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

    // Make Algorithm Object
    GWRGeneralized algorithm;
    algorithm.setCoords(coords);
    algorithm.setDependentVariable(y);
    algorithm.setIndependentVariables(x);
    algorithm.setSpatialWeight(spatial);
    algorithm.setHasHatMatrix(hatmatrix);
    algorithm.setHasIntercept(intercept);
    algorithm.setIsAutoselectBandwidth(optim_bw);
    algorithm.setBandwidthSelectionCriterionType(GWRGeneralized::BandwidthSelectionCriterionType(optim_bw_criterion));

    algorithm.setFamily(GWRGeneralized::Family(size_t(family)));

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
    vec yhat = arma::sum(betas % x, 1);

    List diagnostic =  List::create(
        Named("RSS") = algorithm.getDiagnostic().RSS,
        Named("AIC") = algorithm.getDiagnostic().AIC,
        Named("AICc") = algorithm.getDiagnostic().AICc,
        Named("RSquare") = algorithm.getDiagnostic().RSquare
    );

    List result_list = List::create(
        Named("betas") = betas,
        Named("yhat") = yhat,
        Named("diagnostic") = diagnostic
    );
    if (optim_bw)
    {
        double bw_value = algorithm.spatialWeight().weight<BandwidthWeight>()->bandwidth();
        result_list["bandwidth"] = wrap(bw_value);
    }

    return result_list;
}


// [[Rcpp::export]]
arma::mat gwr_generalized_predict(
    const arma::mat& pcoords,
    const arma::mat& x,
    const arma::vec& y,
    const arma::mat& coords,
    size_t family,
    double bw,
    bool adaptive,
    size_t kernel,
    bool longlat,
    double p,
    double theta,
    bool hatmatrix,
    bool intercept,
    size_t parallel_type,
    const IntegerVector& parallel_arg
) {
    // Convert data types
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

    // Make Algorithm Object
    GWRGeneralized algorithm;
    algorithm.setCoords(coords);
    algorithm.setDependentVariable(y);
    algorithm.setIndependentVariables(x);
    algorithm.setSpatialWeight(spatial);
    algorithm.setHasHatMatrix(hatmatrix);
    algorithm.setHasIntercept(intercept);

    algorithm.setFamily(GWRGeneralized::Family(size_t(family)));

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

    // Return Results
    mat betas;
    try
    {
        betas = algorithm.predict(pcoords);
    }
    catch(const std::exception& e)
    {
        stop(e.what());
    }

    return betas;
}