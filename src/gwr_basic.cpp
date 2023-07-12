#include <Rcpp.h>
#include <armadillo>
#include "utils.h"
#include "gwmodel.h"

using namespace std;
using namespace Rcpp;
using namespace arma;
using namespace gwm;

// [[Rcpp::export]]
List gwr_basic_fit(
    const NumericMatrix& x, const NumericVector& y, const NumericMatrix& coords,
    double bw, bool adaptive, size_t kernel, 
    bool longlat, double p, double theta,
    bool hatmatrix, bool intercept,
    size_t parallel_type, const IntegerVector& parallel_arg,
    bool optim_bw, size_t optim_bw_criterion,
    bool select_model, size_t select_model_criterion, size_t select_model_threshold,
    const CharacterVector& variable_names, bool verbose
) {
    // Convert data types
    mat mx = myas(x);
    vec my = myas(y);
    mat mcoords = myas(coords);
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
    GWRBasic algorithm(mx, my, mcoords, spatial, hatmatrix, intercept);
    algorithm.setIsAutoselectIndepVars(select_model);
    algorithm.setIndepVarSelectionThreshold(select_model_threshold);
    algorithm.setIsAutoselectBandwidth(optim_bw);
    algorithm.setBandwidthSelectionCriterion(GWRBasic::BandwidthSelectionCriterionType(size_t(optim_bw_criterion)));
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
    if (verbose)
        algorithm.setTelegram(make_unique<GWRBasicTelegram>(algorithm, as<vector<string>>(variable_names)));
    algorithm.fit();

    // Return Results
    mat betas = algorithm.betas();
    List result_list = List::create(
        Named("betas") = mywrap(betas),
        Named("betasSE") = mywrap(algorithm.betasSE()),
        Named("sTrace") = mywrap(algorithm.sHat()),
        Named("sHat") = mywrap(algorithm.s()),
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
        mat x = mx.cols(VariableForwardSelector::index2uvec(sel_vars, intercept));
        result_list["fitted"] = mywrap(GWRBasic::Fitted(x, betas));
    }
    else
    {
        result_list["fitted"] = mywrap(GWRBasic::Fitted(mx, betas));
    }

    return result_list;
}


// [[Rcpp::export]]
NumericMatrix gwr_basic_predict(
    const NumericMatrix& pcoords,
    const NumericMatrix& x,
    const NumericVector& y,
    const NumericMatrix& coords,
    double bw,
    bool adaptive,
    size_t kernel,
    bool longlat,
    double p,
    double theta,
    bool intercept,
    size_t parallel_type,
    const IntegerVector& parallel_arg,
    bool verbose
) {
    // Convert data types
    mat mpcoords = myas(pcoords);
    mat mx = myas(x);
    vec my = myas(y);
    mat mcoords = myas(coords);
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
    GWRBasic algorithm(mx, my, mcoords, spatial, false, intercept);
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
    if (verbose)
    {
        algorithm.setTelegram(make_unique<GWRBasicTelegram>(algorithm));
    }

    // Return Results
    mat betas = algorithm.predict(mpcoords);

    return mywrap(betas);
}