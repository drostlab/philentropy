philentropy
===========

## Similarity and Distance Quantification between Probability Functions

The `philentropy` package implements fundamental distance and similarity measures to quantify distances between probability density functions.


## Installation

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

## Tutorials 

 - [Introduction](https://github.com/HajkD/philentropy/blob/master/vignettes/Introduction.Rmd)
 - [Distances and Similarity Measures](https://github.com/HajkD/philentropy/blob/master/vignettes/Distances.Rmd)
 - [Information Theory]()
 
## Important Functions

### Distance Measures
* `distance()` : implements 46 fundamental probability distance (or similarity) measures

### Information Theory

* `H()` : Shannon's Entropy H(X)
* `JE()` : Joint-Entropy H(X,Y)
* `CE()` : Conditional-Entropy H(X | Y)
* `MI()` : Shannon's Mutual Information I(X,Y)
* `KL()` : Kullback–Leibler Divergence
* `JSD()` : Jensen-Shannon Divergence
* `gJSD()` : Generalized Jensen-Shannon Divergence

### Correlation Analyses

* `lin.cor()` : Computes linear correlations 

## Discussions and Bug Reports

I would be very happy to learn more about potential improvements of the concepts and functions
provided in this package.

Furthermore, in case you find some bugs or need additional (more flexible) functionality of parts
of this package, please let me know:

https://github.com/HajkD/philentropy/issues

or find me on [twitter: HajkDrost](https://twitter.com/hajkdrost) 



