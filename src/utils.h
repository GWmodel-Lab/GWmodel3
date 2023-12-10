// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "gwmodel.h"
#include "gwmodelpp/Logger.h"

Rcpp::List mywrap(const gwm::RegressionDiagnostic& diagnostic);
Rcpp::List mywrap(const gwm::VariablesCriterionList& criterion_list);

#define MYMPI_COMM_INFO_DECL \
    int iProcess = 0, nProcess = 1;

#define MYMPI_COMM_INFO_GET \
    MPI_Comm_rank(MPI_COMM_WORLD, &iProcess); \
    MPI_Comm_size(MPI_COMM_WORLD, &nProcess);

#define MYMPI_MASTER_BEGIN \
    if (iProcess == 0) {

#define MYMPI_MASTER_END }

#define MYMPI_WORKER_BEGIN \
    if (nProcess == 0) {

#define MYMPI_WORKER_END }
