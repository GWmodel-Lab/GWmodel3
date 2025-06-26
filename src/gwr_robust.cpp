// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "utils.h"
#include "gwmodel.h"
// #include "telegrams/GWRRobustTelegram.h"

using namespace std;
using namespace Rcpp;
using namespace arma;
using namespace gwm;

// [[Rcpp::export]]
List gwr_robust_fit(
    const arma::mat& x,
    const arma::vec& y,
    const arma::mat& coords,
    double bw,
    bool adaptive,
    size_t kernel,
    bool longlat,
    double p,
    double theta,
    double optim_bw_lower,
    double optim_bw_upper,
    bool hatmatrix,
    bool intercept,
    bool filter,
    size_t parallel_type,
    const IntegerVector& parallel_arg,
    bool optim_bw,
    size_t optim_bw_criterion,
    bool select_model,
    size_t select_model_criterion,
    size_t select_model_threshold,
    const CharacterVector& variable_names,
int verbose
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
    GWRRobust algorithm;
    algorithm.setCoords(coords);
    algorithm.setDependentVariable(y);
    algorithm.setIndependentVariables(x);
    algorithm.setSpatialWeight(spatial);
    algorithm.setHasHatMatrix(hatmatrix);
    algorithm.setHasIntercept(intercept);
    algorithm.setFiltered(filter);
    algorithm.setIsAutoselectIndepVars(select_model);
    algorithm.setIndepVarSelectionThreshold(select_model_threshold);
    algorithm.setIsAutoselectBandwidth(optim_bw);
    algorithm.setBandwidthSelectionCriterion(GWRRobust::BandwidthSelectionCriterionType(size_t(optim_bw_criterion)));
    if (optim_bw_lower > 0.0)
        algorithm.setGoldenLowerBounds(optim_bw_lower);
    if (optim_bw_upper < R_PosInf)
        algorithm.setGoldenUpperBounds(optim_bw_upper);
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
    // if (verbose)
    //     algorithm.setTelegram(make_unique<GWRRobustTelegram>(algorithm, as<vector<string>>(variable_names), verbose));
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
    List result_list = List::create(
        Named("betas") = betas,
        Named("betasSE") = algorithm.betasSE(),
        Named("sTrace") = algorithm.sHat(),
        Named("sHat") = algorithm.s(),
        Named("diagnostic") = mywrap(algorithm.diagnostic())
    );
    if (optim_bw)
    {
        double bw_value = algorithm.spatialWeight().weight<BandwidthWeight>()->bandwidth();
        result_list["bandwidth"] = wrap(bw_value);
    }
    if (select_model)
    {
        vector<size_t> sel_vars = algorithm.selectedVariables();
        result_list["variables"] = wrap(sel_vars);
        result_list["model_sel_criterions"] = mywrap(algorithm.indepVarsSelectionCriterionList());
        mat xs = x.cols(VariableForwardSelector::index2uvec(sel_vars, intercept));
        result_list["fitted"] = GWRRobust::Fitted(xs, betas);
    }
    else
    {
        result_list["fitted"] = GWRRobust::Fitted(x, betas);
    }

    return result_list;
}


// [[Rcpp::export]]
arma::mat gwr_robust_predict(
    const arma::mat& pcoords,
    const arma::mat& x,
    const arma::vec& y,
    const arma::mat& coords,
    double bw,
    bool adaptive,
    size_t kernel,
    bool longlat,
    double p,
    double theta,
    bool intercept,
    size_t parallel_type,
    const IntegerVector& parallel_arg,
    int verbose
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
    GWRRobust algorithm;
    algorithm.setCoords(coords);
    algorithm.setDependentVariable(y);
    algorithm.setIndependentVariables(x);
    algorithm.setSpatialWeight(spatial);
    algorithm.setHasHatMatrix(true);
    algorithm.setHasIntercept(true);
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
    // if (verbose)
    // {
    //     algorithm.setTelegram(make_unique<GWRRobustTelegram>(algorithm, verbose));
    // }

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