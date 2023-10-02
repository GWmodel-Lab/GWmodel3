#include <Rcpp.h>
#include <armadillo>
#include "utils.h"
#include "gwmodel.h"
#include "telegrams/GWRBasicTelegram.h"

#ifdef ENABLE_CUDA_SHARED
#include "gwmodelcuda/IGWRBasicGpuTask.h"
#endif // ENABLE_CUDA_SHARED

using namespace std;
using namespace Rcpp;
using namespace arma;
using namespace gwm;

List gwr_basic_fit_cuda(
    const NumericMatrix& x, const NumericVector& y, const NumericMatrix& coords,
    double bw, bool adaptive, size_t kernel, 
    bool longlat, double p, double theta,
    bool hatmatrix, bool intercept,
    int gpuID, int groupSize,
    bool optim_bw, size_t optim_bw_criterion,
    bool select_model, size_t select_model_criterion, size_t select_model_threshold,
    const CharacterVector& variable_names, int verbose
);

NumericMatrix gwr_basic_predict_cuda(
    const NumericMatrix& pcoords,
    const NumericMatrix& x, const NumericVector& y, const NumericMatrix& coords,
    double bw, bool adaptive, size_t kernel, 
    bool longlat, double p, double theta,
    bool intercept,
    int gpuID, int groupSize,
    int verbose
);

// [[Rcpp::export]]
List gwr_basic_fit(
    const NumericMatrix& x,
    const NumericVector& y,
    const NumericMatrix& coords,
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

#ifdef ENABLE_CUDA_SHARED
    // [For Windows]
    // If parallel type is CUDA, redirect to the specific function
    if (ParallelType(size_t(parallel_type)) == ParallelType::CUDA)
    {
        if (vpar_args.size() < 2) throw std::length_error("CUDA parallelisation needs two parallel args.");
        return gwr_basic_fit_cuda(
            x, y, coords,
            bw, adaptive, kernel,
            longlat, p, theta,
            hatmatrix, intercept,
            vpar_args[0], vpar_args[1],
            optim_bw, optim_bw_criterion,
            select_model, select_model_criterion, select_model_threshold,
            variable_names, verbose
        );
    }
#endif // ENABLE_CUDA_SHARED

    // Convert data types
    mat mx = myas(x);
    vec my = myas(y);
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

    // Make Algorithm Object
    GWRBasic algorithm(mx, my, mcoords, spatial, hatmatrix, intercept);
    algorithm.setIsAutoselectIndepVars(select_model);
    algorithm.setIndepVarSelectionThreshold(select_model_threshold);
    algorithm.setIsAutoselectBandwidth(optim_bw);
    algorithm.setBandwidthSelectionCriterion(GWRBasic::BandwidthSelectionCriterionType(size_t(optim_bw_criterion)));
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
#ifdef ENABLE_CUDA
    case ParallelType::CUDA:
        if (vpar_args.size() < 2) throw std::length_error("CUDA parallelisation needs two parallel args.");
        algorithm.setPrallelType(ParallelType::CUDA);
        algorithm.setGPUId(vpar_args[0]);
        algorithm.setGroupSize(vpar_args[1]);
#endif // ENABLE_CUDA
    default:
        algorithm.setParallelType(ParallelType::SerialOnly);
        break;
    }
    if (verbose)
        algorithm.setTelegram(make_unique<GWRBasicTelegram>(algorithm, as<vector<string>>(variable_names), verbose));
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
    int verbose
) {
    // Convert data types
    vector<int> vpar_args = as< vector<int> >(IntegerVector(parallel_arg));

#ifdef ENABLE_CUDA_SHARED
    // [For Windows]
    // If parallel type is CUDA, redirect to the specific function
    if (ParallelType(size_t(parallel_type)) == ParallelType::CUDA)
    {
        if (vpar_args.size() < 2) throw std::length_error("CUDA parallelisation needs two parallel args.");
        return gwr_basic_predict_cuda(
            pcoords,
            x, y, coords,
            bw, adaptive, kernel,
            longlat, p, theta,
            intercept,
            vpar_args[0], vpar_args[1],
            verbose
        );
    }
#endif // ENABLE_CUDA_SHARED

    mat mpcoords = myas(pcoords);
    mat mx = myas(x);
    vec my = myas(y);
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
#ifdef ENABLE_CUDA
    case ParallelType::CUDA:
        if (vpar_args.size() < 2) throw std::length_error("CUDA parallelisation needs two parallel args.");
        algorithm.setPrallelType(ParallelType::CUDA);
        algorithm.setGPUId(vpar_args[0]);
        algorithm.setGroupSize(vpar_args[1]);
#endif // ENABLE_CUDA
    default:
        algorithm.setParallelType(ParallelType::SerialOnly);
        break;
    }
    if (verbose)
    {
        algorithm.setTelegram(make_unique<GWRBasicTelegram>(algorithm, verbose));
    }

    // Return Results
    mat betas;
    try
    {
        betas = algorithm.predict(mpcoords);
    }
    catch(const std::exception& e)
    {
        stop(e.what());
    }

    return mywrap(betas);
}

List gwr_basic_fit_cuda(
    const NumericMatrix& x, const NumericVector& y, const NumericMatrix& coords,
    double bw, bool adaptive, size_t kernel, 
    bool longlat, double p, double theta,
    bool hatmatrix, bool intercept,
    int gpuID, int groupSize,
    bool optim_bw, size_t optim_bw_criterion,
    bool select_model, size_t select_model_criterion, size_t select_model_threshold,
    const CharacterVector& variable_names, int verbose
) {
#ifdef ENABLE_CUDA_SHARED
    mat mx = myas(x);
    vec my = myas(y);
    mat mcoords = myas(coords);
    
    int distanceType = 0;
    if (longlat)
    {
        distanceType = gwm::Distance::DistanceType::CRSDistance;
    }
    else
    {
        if (p == 2.0 && theta == 0.0)
        {
            distanceType = gwm::Distance::DistanceType::CRSDistance;
        }
        else
        {
            distanceType = gwm::Distance::DistanceType::MinkwoskiDistance;
        }
    }
    if (verbose > 0) Rcout << "** CUDA create task ...";
    auto algorithm = GWRBasicGpuTaskFit_Create(mcoords.n_rows, mx.n_cols, distanceType);
    if (verbose > 0) Rcout << " done\n";
    // Set data
    if (verbose > 0) Rcout << "** CUDA set data ...";
    for (size_t j = 0; j < mx.n_cols; j++)
    {
        for (size_t i = 0; i < mx.n_rows; i++)
        {
            algorithm->setX(i, j, mx(i, j));
        }
    }
    for (size_t i = 0; i < my.n_rows; i++)
    {
        algorithm->setY(i, my(i));
    }
    for (size_t j = 0; j < mcoords.n_cols; j++)
    {
        for (size_t i = 0; i < mcoords.n_rows; i++)
        {
            algorithm->setCoords(i, j, mcoords(i, j));
        }
    }
    if (verbose > 0) Rcout << " done\n";
    // Set distance and weights
    switch (distanceType)
    {
    case gwm::Distance::DistanceType::CRSDistance:
        algorithm->setCRSDistanceGergraphic(longlat);
        break;
    case gwm::Distance::DistanceType::MinkwoskiDistance:
        algorithm->setMinkwoskiDistancePoly(p);
        algorithm->setMinkwoskiDistanceTheta(theta);
        break;
    default:
        algorithm->setCRSDistanceGergraphic(false);
        break;
    }
    algorithm->setBandwidthSize(bw);
    algorithm->setBandwidthAdaptive(adaptive);
    algorithm->setBandwidthKernel(kernel);

    if (optim_bw) {
        algorithm->enableBandwidthOptimization(optim_bw_criterion);
    }
    
    if (select_model) {
        algorithm->enableVariablesOptimization(select_model_threshold);
    }

    if (verbose > 0) Rcout << "** CUDA fit ...";
    if (!algorithm->fit(intercept)) {
        throw std::runtime_error("CUDA did not work successfully.");
    }
    if (verbose > 0) Rcout << " done\n";

    // Get data
    size_t sRows = algorithm->sRows();
    mat betas(size(mx)), betasSE(size(mx)), sHat(sRows, mx.n_rows);
    vec sTrace(2);
    if (verbose > 0) Rcout << "** CUDA get beta ...";
    for (size_t j = 0; j < mx.n_cols; j++)
    {
        for (size_t i = 0; i < mx.n_rows; i++)
        {
            betas(i, j) = algorithm->betas(i, j);
            betasSE(i, j) = algorithm->betasSE(i, j);
        }
    }
    if (verbose > 0) Rcout << " done\n";
    sTrace(0) = algorithm->shat1();
    sTrace(1) = algorithm->shat2();
    for (size_t j = 0; j < mx.n_rows; j++)
    {
        for (size_t i = 0; i < sRows; i++)
        {
            sHat(i, j) = algorithm->s(i, j);
        }   
    }

    RegressionDiagnostic diagnostic;
    diagnostic.RSS = algorithm->diagnosticRSS();
    diagnostic.AIC = algorithm->diagnosticAIC();
    diagnostic.AICc = algorithm->diagnosticAICc();
    diagnostic.ENP = algorithm->diagnosticENP();
    diagnostic.EDF = algorithm->diagnosticEDF();
    diagnostic.RSquare = algorithm->diagnosticRSquare();
    diagnostic.RSquareAdjust = algorithm->diagnosticRSquareAdjust();
    
    // Make result
    List result_list = List::create(
        Named("betas") = mywrap(betas),
        Named("betasSE") = mywrap(betasSE),
        Named("sTrace") = mywrap(sTrace),
        Named("sHat") = mywrap(sHat),
        Named("diagnostic") = mywrap(diagnostic)
    );

    if (optim_bw)
    {
        if (verbose > 0) Rcout << "** CUDA get bw ...";
        result_list["bandwidth"] = wrap(algorithm->optimizedBandwidth());
        if (verbose > 0) Rcout << " done\n";
    }

    if (select_model)
    {
        if (verbose > 0) Rcout << "** CUDA get select model ...";
        vector<size_t> sel_vars(algorithm->selectedVarSize());
        for (size_t i = 0; i < sel_vars.size(); i++)
        {
            sel_vars[i] = algorithm->selectedVar(i);
        }
        result_list["variables"] = wrap(sel_vars);
        
        VariablesCriterionList criterions(algorithm->variableSelectionCriterionSize());
        for (size_t i = 0; i < criterions.size(); i++)
        {
            vector<size_t> vars(algorithm->variableSelectionCriterionItemVarSize(i));
            for (size_t j = 0; j < vars.size(); j++)
            {
                vars[j] = algorithm->variableSelectionCriterionItemVar(i, j);
            }
            criterions[i].first = vars;
            criterions[i].second = algorithm->variableSelectionCriterionItemValue(i);
        }
        result_list["model_sel_criterions"] = mywrap(criterions);
        mat sx = mx.cols(VariableForwardSelector::index2uvec(sel_vars, intercept));
        result_list["fitted"] = mywrap(GWRBasic::Fitted(sx, betas));
        if (verbose > 0) Rcout << " done\n";
    }
    else
    {
        result_list["fitted"] = mywrap(GWRBasic::Fitted(mx, betas));
    }

    if (verbose > 0) Rcout << "** CUDA delete ...";
    GWRBasicGpuTask_Del(algorithm);
    if (verbose > 0) Rcout << " done\n";

    return result_list;
    
#else
    throw std::logic_error("Not supported.");
#endif // ENABLE_CUDA_SHARED
}

NumericMatrix gwr_basic_predict_cuda(
    const NumericMatrix& pcoords,
    const NumericMatrix& x, const NumericVector& y, const NumericMatrix& coords,
    double bw, bool adaptive, size_t kernel, 
    bool longlat, double p, double theta,
    bool intercept,
    int gpuID, int groupSize,
    int verbose
) {
#ifdef ENABLE_CUDA_SHARED
    mat mpcoords = myas(pcoords);
    mat mx = myas(x);
    vec my = myas(y);
    mat mcoords = myas(coords);
    
    int distanceType = 0;
    if (longlat)
    {
        distanceType = gwm::Distance::DistanceType::CRSDistance;
    }
    else
    {
        if (p == 2.0 && theta == 0.0)
        {
            distanceType = gwm::Distance::DistanceType::CRSDistance;
        }
        else
        {
            distanceType = gwm::Distance::DistanceType::MinkwoskiDistance;
        }
    }
    if (verbose > 0) Rcout << "** CUDA create task ...";
    auto algorithm = GWRBasicGpuTaskPredict_Create(mcoords.n_rows, mx.n_cols, distanceType, mpcoords.n_rows);
    if (verbose > 0) Rcout << " done\n";

    // Set data
    if (verbose > 0) Rcout << "** CUDA set data ...";
    for (size_t j = 0; j < mx.n_cols; j++)
    {
        for (size_t i = 0; i < mx.n_rows; i++)
        {
            algorithm->setX(i, j, mx(i, j));
        }
    }
    for (size_t i = 0; i < my.n_rows; i++)
    {
        algorithm->setY(i, my(i));
    }
    for (size_t j = 0; j < mcoords.n_cols; j++)
    {
        for (size_t i = 0; i < mcoords.n_rows; i++)
        {
            algorithm->setCoords(i, j, mcoords(i, j));
        }
    }
    for (size_t j = 0; j < mpcoords.n_cols; j++)
    {
        for (size_t i = 0; i < mpcoords.n_rows; i++)
        {
            algorithm->setPredictLocations(i, j, mpcoords(i, j));
        }
    }
    
    if (verbose > 0) Rcout << " done\n";
    // Set distance and weights
    switch (distanceType)
    {
    case gwm::Distance::DistanceType::CRSDistance:
        algorithm->setCRSDistanceGergraphic(longlat);
        break;
    case gwm::Distance::DistanceType::MinkwoskiDistance:
        algorithm->setMinkwoskiDistancePoly(p);
        algorithm->setMinkwoskiDistanceTheta(theta);
        break;
    default:
        algorithm->setCRSDistanceGergraphic(false);
        break;
    }
    algorithm->setBandwidthSize(bw);
    algorithm->setBandwidthAdaptive(adaptive);
    algorithm->setBandwidthKernel(kernel);

    if (verbose > 0) Rcout << "** CUDA predict ...";
    if (!algorithm->predict(intercept)) {
        throw std::runtime_error("CUDA did not work successfully.");
    }
    if (verbose > 0) Rcout << " done\n";

    // Get data
    mat betas(size(mx));
    if (verbose > 0) Rcout << "** CUDA get beta ...";
    for (size_t j = 0; j < mx.n_cols; j++)
    {
        for (size_t i = 0; i < mx.n_rows; i++)
        {
            betas(i, j) = algorithm->betas(i, j);
        }
    }
    if (verbose > 0) Rcout << " done\n";

    if (verbose > 0) Rcout << "** CUDA delete ...";
    GWRBasicGpuTask_Del(algorithm);
    if (verbose > 0) Rcout << " done\n";

    return mywrap(betas);
    
#else
    throw std::logic_error("Not supported.");
#endif // ENABLE_CUDA_SHARED
}
