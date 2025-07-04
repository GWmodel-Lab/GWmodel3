% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gwr_localcollinearity.R
\name{gwr_lcr}
\alias{gwr_lcr}
\alias{plot.gwlcrm}
\alias{coef.gwlcrm}
\alias{fitted.gwlcrm}
\alias{residuals.gwlcrm}
\title{Calibrate a GWR local collinearity model}
\usage{
gwr_lcr(
  formula,
  data,
  bw = NA,
  adaptive = FALSE,
  kernel = c("gaussian", "exp", "bisquare", "tricube", "boxcar"),
  longlat = FALSE,
  p = 2,
  theta = 0,
  lambda = 0,
  lambda_adjust = FALSE,
  cn_thresh = 30,
  hatmatrix = TRUE,
  parallel_method = c("no", "omp"),
  parallel_arg = c(0),
  verbose = FALSE
)

\method{plot}{gwlcrm}(x, y, ..., columns)

\method{coef}{gwlcrm}(object, ...)

\method{fitted}{gwlcrm}(object, ...)

\method{residuals}{gwlcrm}(object, ...)
}
\arguments{
\item{formula}{Regresison model.}

\item{data}{A \code{sf} objects.}

\item{bw}{Either a value to set the size of bandwidth,
or \code{CV} to set the criterion for bandwidth auto-optimization process.
Note that if \code{NA} or other non-numeric value is setted,
this parameter will be reset to \code{Inf}.}

\item{adaptive}{Whether the bandwidth value is adaptive or not.}

\item{kernel}{Kernel function used.}

\item{longlat}{Whether the coordinates}

\item{p}{Power of the Minkowski distance,
default to 2, i.e., Euclidean distance.}

\item{theta}{Angle in radian to roate the coordinate system, default to 0.}

\item{lambda}{Option for a globally-defined (constant) ridge parameter.
Default is lambda=0, which gives a basic GWR fit}

\item{lambda_adjust}{A locally-varying ridge parameter.Default FALSE, refers to:
-i a basic GWR without a local ridge adjustment
(i.e. lambda=0, everywhere);
-ii a penalised GWR with a global ridge adjustment
(i.e. lambda is user-specified as some constant, other than 0 everywhere);
if TRUE, use cn.tresh to set the maximum condition number.
Here for locations with a condition number (for its local design matrix)
above this user-specified threshold, a local ridge parameter is found}

\item{cn_thresh}{maximum value for condition number, commonly set between 20 and 30}

\item{hatmatrix}{If TRUE, great circle will be caculated.}

\item{parallel_method}{Parallel method, multithreading (\code{omp}) is available}

\item{parallel_arg}{Parallel method argument.}

\item{verbose}{Whether to print additional information.}

\item{x}{A "gwlcrm" object.}

\item{y}{Ignored.}

\item{\dots}{Additional arguments passing to \code{\link[=residuals]{residuals()}}.}

\item{columns}{Column names to plot.
If it is missing or non-character value, all coefficient columns are plottd.}

\item{object}{A "gwlcrm" object.}
}
\value{
A \code{gwlcrm} object.
}
\description{
Calibrate a GWR local collinearity model
}
\section{Functions}{
\itemize{
\item \code{plot(gwlcrm)}: Plot the result of GWR local collinearity model.

\item \code{coef(gwlcrm)}: Get coefficients of a GWR local collinearity model.

\item \code{fitted(gwlcrm)}: Get fitted values of a GWR local collinearity model.

\item \code{residuals(gwlcrm)}: Get residuals of a GWR local collinearity model.

}}
\examples{
data(LondonHP)

# Basic usage
gwr_lcr(PURCHASE ~ FLOORSZ + UNEMPLOY, LondonHP, 64, TRUE)

# Bandwidth Optimization
m <- gwr_lcr(PURCHASE ~ FLOORSZ + UNEMPLOY + PROF, LondonHP, 'CV', TRUE)
m

plot(m)

coef(m)

fitted(m)

residuals(m)

}
