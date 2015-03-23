philentropy
===========

## An Information Theory And Statistics Framework for Predictive Analytics

The `philentropy` Package Provides an Information Theory and Statistics Framework to Perform Predictive Analytics.


## Fast Installation Guide

```r
# install.packages("devtools")

# install the current version of philentropy on your system
library(devtools)
install_github("HajkD/philentropy", build_vignettes = TRUE, dependencies = TRUE)

# On Windows, this won't work - see ?build_github_devtools
install_github("HajkD/philentropy", build_vignettes = TRUE, dependencies = TRUE)

# When working with Windows, first you need to install the
# R package: rtools -> install.packages("rtools")

# Afterwards you can install devtools -> install.packages("devtools")
# and then you can run:

devtools::install_github("HajkD/philentropy", build_vignettes = TRUE, dependencies = TRUE)

# and then call it from the library
library("philentropy", lib.loc = "C:/Program Files/R/R-3.1.1/library")

```

## Discussions and Bug Reports

I would be very happy to learn more about potential improvements of the concepts and functions
provided in this package.

Furthermore, in case you find some bugs or need additional (more flexible) functionality of parts
of this package, please let me know:

hajk-georg.drost@informatik.uni-halle.de




