#include <Rcpp.h>
#include <armadillo>
#include "utils.h"
#include "gwmodel.h"
#include "telegrams/GWDRTelegram.h"

using namespace std;
using namespace Rcpp;
using namespace arma;
using namespace gwm;

// [[Rcpp::export]]
List gwdr_fit(
    const NumericMatrix& x,
    const NumericVector& y,
    const NumericMatrix& coords,
    const NumericVector& bw,
    const LogicalVector& adaptive,
    const IntegerVector& kernel,
    bool intercept,
    bool hatmatrix,
    size_t parallel_type,
    const IntegerVector& parallel_arg,
    bool optim_bw,
    size_t optim_bw_criterion,
    double optim_threashold,
    double optim_step,
    size_t optim_max_iter,
    bool select_model,
    size_t select_model_threshold,
    const CharacterVector& variable_names,
    int verbose
) {
    // Convert data types
    arma::mat mx = myas(x);
    arma::vec my = myas(y);
    arma::mat mcoords = myas(coords);
    std::vector<int> vpar_args = as< std::vector<int> >(Rcpp::IntegerVector(parallel_arg));

    // Make Spatial Weight
    size_t nDim = (size_t)mcoords.n_cols;
    auto vbw = as< vector<double> >(NumericVector(bw));
    auto vadaptive = as< vector<bool> >(LogicalVector(adaptive));
    auto vkernel = as< vector<int> >(IntegerVector(kernel));
    vector<SpatialWeight> spatials;
    for (size_t i = 0; i < nDim; i++)
    {
        BandwidthWeight bandwidth(vbw[i] * mcoords.n_rows, vadaptive[i], BandwidthWeight::KernelFunctionType(vkernel[i]));
        OneDimDistance distance;
        spatials.push_back(SpatialWeight(&bandwidth, &distance));
    }
    
    // Make Algorithm Object
    GWDR algorithm(mx, my, mcoords, spatials);
    algorithm.setHasHatMatrix(hatmatrix);
    algorithm.setBandwidthCriterionType(GWDR::BandwidthCriterionType(size_t(optim_bw_criterion)));

    if (select_model)
    {
        algorithm.setEnableIndepVarSelect(true);
        algorithm.setIndepVarSelectThreshold(select_model_threshold);
    }

    if (optim_bw) 
    {
        algorithm.setEnableBandwidthOptimize(true);
        algorithm.setBandwidthOptimizeEps(optim_threashold);
        algorithm.setBandwidthOptimizeStep(optim_step);
        algorithm.setBandwidthOptimizeMaxIter(optim_max_iter);
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
    if (verbose > 0)
    {
        algorithm.setTelegram(make_unique<GWDRTelegram>(algorithm, as<vector<string>>(variable_names), verbose));
    }
    algorithm.fit();
    
    // Return Results
    mat betas = algorithm.betas();
    vec fitted = sum(mx % betas, 1);
    List result_list = List::create(
        Named("betas") = mywrap(betas),
        Named("betasSE") = mywrap(algorithm.betasSE()),
        Named("sTrace") = mywrap(algorithm.sHat()),
        Named("sHat") = mywrap(algorithm.s()),
        Named("diagnostic") = mywrap(algorithm.diagnostic())
    );
    if (optim_bw)
    {
        vector<double> bw_value;
        const vector<SpatialWeight>& spatialWeights = algorithm.spatialWeights();
        for (size_t i = 0; i < nDim; i++)
        {
            bw_value.push_back(spatialWeights[i].weight<BandwidthWeight>()->bandwidth() / double(mcoords.n_rows));
        }
        result_list["bw_value"] = wrap(bw_value);
    }
    if (select_model)
    {
        vector<size_t> sel_vars = algorithm.selectedVariables();
        result_list["variables"] = wrap(sel_vars);
        result_list["model_sel_criterions"] = mywrap(algorithm.indepVarCriterionList());
        mat x = mx.cols(VariableForwardSelector::index2uvec(sel_vars, intercept));
        result_list["fitted"] = mywrap(GWRBasic::Fitted(x, betas));
    }
    else
    {
        result_list["fitted"] = mywrap(GWRBasic::Fitted(mx, betas));
    }

    return result_list;
}
