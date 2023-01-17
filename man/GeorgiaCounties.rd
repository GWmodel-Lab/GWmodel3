\name{GeorgiaCounties}
\alias{GeorgiaCounties}
\docType{data}
\title{Georgia counties data (sf)}
\description{
  The Georgia census data with boundaries for mapping
}
\usage{data(GeorgiaCounties)}
\details{
This data set can also be found in GWR 3 and in spgwr.
}
\examples{
library(sf)
data(GeorgiaCensus)
data(GeorgiaCounties)
georgia_sf <- st_as_sf(GeorgiaCensus, coords = c("X", "Y"))
plot(GeorgiaCounties$geometry)
plot(georgia_sf["TotPop90"], add = TRUE)
}
\keyword{data}
\concept{Georgia counties}
