// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// gwdr_fit
List gwdr_fit(const arma::mat& x, const arma::vec& y, const arma::mat& coords, const NumericVector& bw, const LogicalVector& adaptive, const IntegerVector& kernel, bool intercept, bool hatmatrix, size_t parallel_type, const IntegerVector& parallel_arg, bool optim_bw, size_t optim_bw_criterion, double optim_threashold, double optim_step, size_t optim_max_iter, bool select_model, size_t select_model_threshold, const CharacterVector& variable_names, int verbose);
RcppExport SEXP _GWmodel3_gwdr_fit(SEXP xSEXP, SEXP ySEXP, SEXP coordsSEXP, SEXP bwSEXP, SEXP adaptiveSEXP, SEXP kernelSEXP, SEXP interceptSEXP, SEXP hatmatrixSEXP, SEXP parallel_typeSEXP, SEXP parallel_argSEXP, SEXP optim_bwSEXP, SEXP optim_bw_criterionSEXP, SEXP optim_threasholdSEXP, SEXP optim_stepSEXP, SEXP optim_max_iterSEXP, SEXP select_modelSEXP, SEXP select_model_thresholdSEXP, SEXP variable_namesSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type coords(coordsSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type bw(bwSEXP);
    Rcpp::traits::input_parameter< const LogicalVector& >::type adaptive(adaptiveSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type kernel(kernelSEXP);
    Rcpp::traits::input_parameter< bool >::type intercept(interceptSEXP);
    Rcpp::traits::input_parameter< bool >::type hatmatrix(hatmatrixSEXP);
    Rcpp::traits::input_parameter< size_t >::type parallel_type(parallel_typeSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type parallel_arg(parallel_argSEXP);
    Rcpp::traits::input_parameter< bool >::type optim_bw(optim_bwSEXP);
    Rcpp::traits::input_parameter< size_t >::type optim_bw_criterion(optim_bw_criterionSEXP);
    Rcpp::traits::input_parameter< double >::type optim_threashold(optim_threasholdSEXP);
    Rcpp::traits::input_parameter< double >::type optim_step(optim_stepSEXP);
    Rcpp::traits::input_parameter< size_t >::type optim_max_iter(optim_max_iterSEXP);
    Rcpp::traits::input_parameter< bool >::type select_model(select_modelSEXP);
    Rcpp::traits::input_parameter< size_t >::type select_model_threshold(select_model_thresholdSEXP);
    Rcpp::traits::input_parameter< const CharacterVector& >::type variable_names(variable_namesSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(gwdr_fit(x, y, coords, bw, adaptive, kernel, intercept, hatmatrix, parallel_type, parallel_arg, optim_bw, optim_bw_criterion, optim_threashold, optim_step, optim_max_iter, select_model, select_model_threshold, variable_names, verbose));
    return rcpp_result_gen;
END_RCPP
}
// gwr_basic_fit
List gwr_basic_fit(const arma::mat& x, const arma::vec& y, const arma::mat& coords, double bw, bool adaptive, size_t kernel, bool longlat, double p, double theta, double optim_bw_lower, double optim_bw_upper, bool hatmatrix, bool intercept, size_t parallel_type, const IntegerVector& parallel_arg, bool optim_bw, size_t optim_bw_criterion, bool select_model, size_t select_model_criterion, size_t select_model_threshold, const CharacterVector& variable_names, int verbose);
RcppExport SEXP _GWmodel3_gwr_basic_fit(SEXP xSEXP, SEXP ySEXP, SEXP coordsSEXP, SEXP bwSEXP, SEXP adaptiveSEXP, SEXP kernelSEXP, SEXP longlatSEXP, SEXP pSEXP, SEXP thetaSEXP, SEXP optim_bw_lowerSEXP, SEXP optim_bw_upperSEXP, SEXP hatmatrixSEXP, SEXP interceptSEXP, SEXP parallel_typeSEXP, SEXP parallel_argSEXP, SEXP optim_bwSEXP, SEXP optim_bw_criterionSEXP, SEXP select_modelSEXP, SEXP select_model_criterionSEXP, SEXP select_model_thresholdSEXP, SEXP variable_namesSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type coords(coordsSEXP);
    Rcpp::traits::input_parameter< double >::type bw(bwSEXP);
    Rcpp::traits::input_parameter< bool >::type adaptive(adaptiveSEXP);
    Rcpp::traits::input_parameter< size_t >::type kernel(kernelSEXP);
    Rcpp::traits::input_parameter< bool >::type longlat(longlatSEXP);
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< double >::type optim_bw_lower(optim_bw_lowerSEXP);
    Rcpp::traits::input_parameter< double >::type optim_bw_upper(optim_bw_upperSEXP);
    Rcpp::traits::input_parameter< bool >::type hatmatrix(hatmatrixSEXP);
    Rcpp::traits::input_parameter< bool >::type intercept(interceptSEXP);
    Rcpp::traits::input_parameter< size_t >::type parallel_type(parallel_typeSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type parallel_arg(parallel_argSEXP);
    Rcpp::traits::input_parameter< bool >::type optim_bw(optim_bwSEXP);
    Rcpp::traits::input_parameter< size_t >::type optim_bw_criterion(optim_bw_criterionSEXP);
    Rcpp::traits::input_parameter< bool >::type select_model(select_modelSEXP);
    Rcpp::traits::input_parameter< size_t >::type select_model_criterion(select_model_criterionSEXP);
    Rcpp::traits::input_parameter< size_t >::type select_model_threshold(select_model_thresholdSEXP);
    Rcpp::traits::input_parameter< const CharacterVector& >::type variable_names(variable_namesSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(gwr_basic_fit(x, y, coords, bw, adaptive, kernel, longlat, p, theta, optim_bw_lower, optim_bw_upper, hatmatrix, intercept, parallel_type, parallel_arg, optim_bw, optim_bw_criterion, select_model, select_model_criterion, select_model_threshold, variable_names, verbose));
    return rcpp_result_gen;
END_RCPP
}
// gwr_basic_predict
arma::mat gwr_basic_predict(const arma::mat& pcoords, const arma::mat& x, const arma::vec& y, const arma::mat& coords, double bw, bool adaptive, size_t kernel, bool longlat, double p, double theta, bool intercept, size_t parallel_type, const IntegerVector& parallel_arg, int verbose);
RcppExport SEXP _GWmodel3_gwr_basic_predict(SEXP pcoordsSEXP, SEXP xSEXP, SEXP ySEXP, SEXP coordsSEXP, SEXP bwSEXP, SEXP adaptiveSEXP, SEXP kernelSEXP, SEXP longlatSEXP, SEXP pSEXP, SEXP thetaSEXP, SEXP interceptSEXP, SEXP parallel_typeSEXP, SEXP parallel_argSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type pcoords(pcoordsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type coords(coordsSEXP);
    Rcpp::traits::input_parameter< double >::type bw(bwSEXP);
    Rcpp::traits::input_parameter< bool >::type adaptive(adaptiveSEXP);
    Rcpp::traits::input_parameter< size_t >::type kernel(kernelSEXP);
    Rcpp::traits::input_parameter< bool >::type longlat(longlatSEXP);
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< bool >::type intercept(interceptSEXP);
    Rcpp::traits::input_parameter< size_t >::type parallel_type(parallel_typeSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type parallel_arg(parallel_argSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(gwr_basic_predict(pcoords, x, y, coords, bw, adaptive, kernel, longlat, p, theta, intercept, parallel_type, parallel_arg, verbose));
    return rcpp_result_gen;
END_RCPP
}
// gwr_multiscale_fit
List gwr_multiscale_fit(const arma::mat& x, const arma::vec& y, const arma::mat& coords, const NumericVector& bw, const LogicalVector& adaptive, const IntegerVector& kernel, const LogicalVector& longlat, const NumericVector& p, const NumericVector& theta, const LogicalVector& optim_bw, const IntegerVector& optim_bw_criterion, const NumericVector& threashold, const IntegerVector& initial_type, const LogicalVector& centered, double optim_bw_lower, double optim_bw_upper, size_t criterion, bool hatmatrix, bool intercept, size_t retry_times, size_t max_iterations, size_t parallel_type, const IntegerVector& parallel_arg, const CharacterVector& variable_names, int verbose);
RcppExport SEXP _GWmodel3_gwr_multiscale_fit(SEXP xSEXP, SEXP ySEXP, SEXP coordsSEXP, SEXP bwSEXP, SEXP adaptiveSEXP, SEXP kernelSEXP, SEXP longlatSEXP, SEXP pSEXP, SEXP thetaSEXP, SEXP optim_bwSEXP, SEXP optim_bw_criterionSEXP, SEXP threasholdSEXP, SEXP initial_typeSEXP, SEXP centeredSEXP, SEXP optim_bw_lowerSEXP, SEXP optim_bw_upperSEXP, SEXP criterionSEXP, SEXP hatmatrixSEXP, SEXP interceptSEXP, SEXP retry_timesSEXP, SEXP max_iterationsSEXP, SEXP parallel_typeSEXP, SEXP parallel_argSEXP, SEXP variable_namesSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type coords(coordsSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type bw(bwSEXP);
    Rcpp::traits::input_parameter< const LogicalVector& >::type adaptive(adaptiveSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type kernel(kernelSEXP);
    Rcpp::traits::input_parameter< const LogicalVector& >::type longlat(longlatSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type p(pSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const LogicalVector& >::type optim_bw(optim_bwSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type optim_bw_criterion(optim_bw_criterionSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type threashold(threasholdSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type initial_type(initial_typeSEXP);
    Rcpp::traits::input_parameter< const LogicalVector& >::type centered(centeredSEXP);
    Rcpp::traits::input_parameter< double >::type optim_bw_lower(optim_bw_lowerSEXP);
    Rcpp::traits::input_parameter< double >::type optim_bw_upper(optim_bw_upperSEXP);
    Rcpp::traits::input_parameter< size_t >::type criterion(criterionSEXP);
    Rcpp::traits::input_parameter< bool >::type hatmatrix(hatmatrixSEXP);
    Rcpp::traits::input_parameter< bool >::type intercept(interceptSEXP);
    Rcpp::traits::input_parameter< size_t >::type retry_times(retry_timesSEXP);
    Rcpp::traits::input_parameter< size_t >::type max_iterations(max_iterationsSEXP);
    Rcpp::traits::input_parameter< size_t >::type parallel_type(parallel_typeSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type parallel_arg(parallel_argSEXP);
    Rcpp::traits::input_parameter< const CharacterVector& >::type variable_names(variable_namesSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(gwr_multiscale_fit(x, y, coords, bw, adaptive, kernel, longlat, p, theta, optim_bw, optim_bw_criterion, threashold, initial_type, centered, optim_bw_lower, optim_bw_upper, criterion, hatmatrix, intercept, retry_times, max_iterations, parallel_type, parallel_arg, variable_names, verbose));
    return rcpp_result_gen;
END_RCPP
}
// gwss_fit
List gwss_fit(const NumericMatrix& x, const NumericMatrix& coords, int mode, bool quantile, double bw, bool adaptive, size_t kernel, bool longlat, double p, double theta, size_t parallel_type, const IntegerVector& parallel_arg);
RcppExport SEXP _GWmodel3_gwss_fit(SEXP xSEXP, SEXP coordsSEXP, SEXP modeSEXP, SEXP quantileSEXP, SEXP bwSEXP, SEXP adaptiveSEXP, SEXP kernelSEXP, SEXP longlatSEXP, SEXP pSEXP, SEXP thetaSEXP, SEXP parallel_typeSEXP, SEXP parallel_argSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type coords(coordsSEXP);
    Rcpp::traits::input_parameter< int >::type mode(modeSEXP);
    Rcpp::traits::input_parameter< bool >::type quantile(quantileSEXP);
    Rcpp::traits::input_parameter< double >::type bw(bwSEXP);
    Rcpp::traits::input_parameter< bool >::type adaptive(adaptiveSEXP);
    Rcpp::traits::input_parameter< size_t >::type kernel(kernelSEXP);
    Rcpp::traits::input_parameter< bool >::type longlat(longlatSEXP);
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< size_t >::type parallel_type(parallel_typeSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type parallel_arg(parallel_argSEXP);
    rcpp_result_gen = Rcpp::wrap(gwss_fit(x, coords, mode, quantile, bw, adaptive, kernel, longlat, p, theta, parallel_type, parallel_arg));
    return rcpp_result_gen;
END_RCPP
}
// gwss_average
List gwss_average(const NumericMatrix& x, const NumericMatrix& coords, bool quantile, double bw, bool adaptive, size_t kernel, bool longlat, double p, double theta, size_t parallel_type, const IntegerVector& parallel_arg);
RcppExport SEXP _GWmodel3_gwss_average(SEXP xSEXP, SEXP coordsSEXP, SEXP quantileSEXP, SEXP bwSEXP, SEXP adaptiveSEXP, SEXP kernelSEXP, SEXP longlatSEXP, SEXP pSEXP, SEXP thetaSEXP, SEXP parallel_typeSEXP, SEXP parallel_argSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type coords(coordsSEXP);
    Rcpp::traits::input_parameter< bool >::type quantile(quantileSEXP);
    Rcpp::traits::input_parameter< double >::type bw(bwSEXP);
    Rcpp::traits::input_parameter< bool >::type adaptive(adaptiveSEXP);
    Rcpp::traits::input_parameter< size_t >::type kernel(kernelSEXP);
    Rcpp::traits::input_parameter< bool >::type longlat(longlatSEXP);
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< size_t >::type parallel_type(parallel_typeSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type parallel_arg(parallel_argSEXP);
    rcpp_result_gen = Rcpp::wrap(gwss_average(x, coords, quantile, bw, adaptive, kernel, longlat, p, theta, parallel_type, parallel_arg));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_GWmodel3_gwdr_fit", (DL_FUNC) &_GWmodel3_gwdr_fit, 19},
    {"_GWmodel3_gwr_basic_fit", (DL_FUNC) &_GWmodel3_gwr_basic_fit, 22},
    {"_GWmodel3_gwr_basic_predict", (DL_FUNC) &_GWmodel3_gwr_basic_predict, 14},
    {"_GWmodel3_gwr_multiscale_fit", (DL_FUNC) &_GWmodel3_gwr_multiscale_fit, 25},
    {"_GWmodel3_gwss_fit", (DL_FUNC) &_GWmodel3_gwss_fit, 12},
    {"_GWmodel3_gwss_average", (DL_FUNC) &_GWmodel3_gwss_average, 11},
    {NULL, NULL, 0}
};

RcppExport void R_init_GWmodel3(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
