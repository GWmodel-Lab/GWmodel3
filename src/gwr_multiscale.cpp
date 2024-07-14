// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "utils.h"
#include "gwmodel.h"
#include "telegrams/GWRMultiscaleTelegram.h"

#ifdef ENABLE_MPI
#include <mpi.h>
#endif // ENABLE_MPI

using namespace std;
using namespace Rcpp;
using namespace arma;
using namespace gwm;

// [[Rcpp::export]]
List gwr_multiscale_fit (
    const arma::mat& x,
    const arma::vec& y,
    const arma::mat& coords,
    const NumericVector& bw,
    const LogicalVector& adaptive,
    const IntegerVector& kernel,
    const LogicalVector& longlat,
    const NumericVector& p,
    const NumericVector& theta,
    const LogicalVector& optim_bw,
    const IntegerVector& optim_bw_criterion,
    const NumericVector& threashold,
    const IntegerVector& initial_type,
    const LogicalVector& centered,
    double optim_bw_lower,
    double optim_bw_upper,
    size_t criterion,
    bool hatmatrix,
    bool intercept,
    size_t retry_times,
    size_t max_iterations,
    size_t parallel_type,
    const IntegerVector& parallel_arg,
    const CharacterVector& variable_names,
    int verbose
) {
    // Make Spatial Weight
    vector<int> vpar_args = as< vector<int> >(IntegerVector(parallel_arg));
    size_t nVar = (size_t)x.n_cols;
    auto vbw = as< vector<double> >(NumericVector(bw));
    auto vadaptive = as< vector<bool> >(LogicalVector(adaptive));
    auto vkernel = as< vector<int> >(IntegerVector(kernel));
    auto vlonglat = as< vector<bool> >(LogicalVector(longlat));
    auto vp = as< vector<double> >(NumericVector(p));
    auto vtheta = as< vector<double> >(NumericVector(theta));
    auto voptim_bw = as< vector<bool> >(LogicalVector(optim_bw));
    auto voptim_bw_criterion = as< vector<int> >(IntegerVector(optim_bw_criterion));
    auto vinitial_type = as< vector<int> >(IntegerVector(initial_type));
    auto vcentered = as< vector<bool> >(LogicalVector(centered));
    auto vthreshold = as< vector<double> >(NumericVector(threashold));
    vector<SpatialWeight> spatials;
    for (size_t i = 0; i < nVar; i++)
    {
        BandwidthWeight bandwidth(vbw[i], vadaptive[i], BandwidthWeight::KernelFunctionType(vkernel[i]));
        Distance* distance = nullptr;
        if (vlonglat[i]) distance = new CRSDistance(true);
        else
        {
            if (vp[i] == 2.0 && vtheta[i] == 0.0) distance = new CRSDistance(false);
            else distance = new MinkwoskiDistance(vp[i], vtheta[i]);
        }
        spatials.push_back(SpatialWeight(&bandwidth, distance));
    }
    vector<GWRMultiscale::BandwidthInitilizeType> bandwidthInitialize(vinitial_type.size());
    transform(vinitial_type.begin(), vinitial_type.end(), bandwidthInitialize.begin(), [](int x) {
        return GWRMultiscale::BandwidthInitilizeType(x);
    });
    vector<GWRMultiscale::BandwidthSelectionCriterionType> bandwidthSelectionApproach(voptim_bw_criterion.size());
    transform(voptim_bw_criterion.begin(), voptim_bw_criterion.end(), bandwidthSelectionApproach.begin(), [](int x) {
        return GWRMultiscale::BandwidthSelectionCriterionType(x);
    });
    
    // Make Algorithm Object
    GWRMultiscale algorithm(x, y, coords, spatials);
    algorithm.setIndependentVariables(x);
    algorithm.setDependentVariable(y);
    algorithm.setCoords(coords);
    algorithm.setSpatialWeights(spatials);
    algorithm.setPreditorCentered(vcentered);
    algorithm.setBandwidthInitilize(bandwidthInitialize);
    algorithm.setBandwidthSelectionApproach(bandwidthSelectionApproach);
    algorithm.setBandwidthSelectThreshold(vthreshold);
    if (optim_bw_lower > 0.0)
        algorithm.setGoldenLowerBounds(optim_bw_lower);
    if (optim_bw_upper < R_PosInf)
        algorithm.setGoldenUpperBounds(optim_bw_upper);
    algorithm.setCriterionType(GWRMultiscale::BackFittingCriterionType(size_t(criterion)));
    algorithm.setHasHatMatrix(hatmatrix);
    algorithm.setBandwidthSelectRetryTimes(retry_times);
    algorithm.setMaxIteration(max_iterations);
    algorithm.setParallelType(ParallelType(parallel_type));
    switch (ParallelType(size_t(parallel_type)))
    {
    case ParallelType::OpenMP:
#ifdef ENABLE_MPI
    case ParallelType::MPI_MP:
#endif // ENABLE_MPI
#ifdef _OPENMP
        algorithm.setParallelType(ParallelType::OpenMP);
        algorithm.setOmpThreadNum(vpar_args[0]);
        break;
#else
        throw std::logic_error("OpenMP method not implemented.");
#endif
    default:
        break;
    }
    
    MYMPI_COMM_INFO_DECL
#ifdef ENABLE_MPI
    if (algorithm.parallelType() & ParallelType::MPI)
    {
        MYMPI_COMM_INFO_GET
        algorithm.setWorkerId(iProcess);
        algorithm.setWorkerNum(nProcess);
        MYMPI_MASTER_BEGIN
        Rcout << "* MPI mode\n";
        MYMPI_MASTER_END
    }
#endif // ENABLE_MPI

    MYMPI_MASTER_BEGIN
    if (verbose > 0)
    {
        algorithm.setTelegram(make_unique<GWRMultiscaleTelegram>(algorithm, as<vector<string>>(variable_names), verbose));
    }
    MYMPI_MASTER_END

    try
    {
        algorithm.fit();
    }
    catch(const std::exception& e)
    {
        stop(e.what());
    }

    List result_list;
    
    MYMPI_MASTER_BEGIN
    // Get bandwidth
    vector<double> bw_value;
    const vector<SpatialWeight>& spatialWeights = algorithm.spatialWeights();
    for (size_t i = 0; i < nVar; i++)
    {
        bw_value.push_back(spatialWeights[i].weight<BandwidthWeight>()->bandwidth());
    }
    
    // Return Results
    mat betas = algorithm.betas();
    vec fitted = sum(x % betas, 1);
    result_list = List::create(
        Named("betas") = betas,
        Named("diagnostic") = mywrap(algorithm.diagnostic()),
        Named("bw_value") = wrap(bw_value),
        Named("fitted") = fitted
    );
    MYMPI_MASTER_END

#ifdef ENABLE_MPI
    if (parallel_type & ParallelType::MPI)
        result_list["mpi_rank"] = iProcess;
#endif // ENABLE_MPI

    return result_list;
}
