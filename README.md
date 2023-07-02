# GWmodel3: the next-generation R package for geographically weighted modeling

## Overview

This package includes techniques from a particular branch of spatial statistics, 
termed geographically weighted (GW) models.
GW models suit situations when data are not described well by some global model,
but where there are spatial regions where a suitably localized calibration provides a better description.

From version 3.0, the goal of **GWmodel** is to provide more conscious and easier user interfaces
and high-performance implementations by refactoring R functions and internal C++ code.
And the package [sf](https://r-spatial.github.io/sf/) now is set as the default dependency to manipulate spatial data.
We believe with the newly designed interfaces and underlying code,
users will get fluent and highly consistent experiences.

## Installation

Install from GitHub:

```R
devtools::install_github("GWmodel-Lab/GWmodel3")
```

## Getting started

```R
library(tidyr)
```

Now the following models have been implemented in this package:

- [x] GW regression of various forms ([Brunsdon et al., 1996](https://doi.org/10.1111/j.1538-4632.1996.tb00936.x))
- [ ] GW summary statistics ([Brunsdon et al., 2002](https://doi.org/10.1016/s0198-9715(01)00009-6))
- [ ] GW principal components analysis ([Harris et al., 2011](https://doi.org/10.1080/13658816.2011.554838))
- [ ] GW discriminant analysis ([Brunsdon et al., 2007](https://doi.org/10.1111/j.1538-4632.2007.00709.x))

Please find vignettes for more information.

## Related work

From version 3.0, **GWmodel** is based on a pure C++ library --- [libgwmodel](https://github.com/GWmodel-Lab/libgwmdoel).
This library implements all models, and **GWmodel** just calls this package by translating inputs and outputs.
Besides, **GWmodel** also provides some handy functions for the convenience of R users.

For the old version, please turn to the [old repository](https://github.com/GWmodel-Lab/GWmodel).

## Getting help

If you encounter a bug, please create an issue [here](https://github.com/GWmodel-Lab/GWmodel3/issues).
It would be better for us if a minimal reproducible example is also provided.
