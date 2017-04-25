philentropy
===========

### Similarity and Distance Quantification between Probability Functions

> Describe and understand the world through data.

Data collection and data comparison are the foundations of scientific research.
_Mathematics_ provides the abstract framework to describe patterns we observe in nature and _Statistics_ provides the
framework to quantify the uncertainty of these patterns. In statistics natural patterns
are described in form of probability distributions which either follow a fixed pattern (parametric distributions) or more dynamic patterns (non-parametric distributions).

The `philentropy` package implements fundamental distance and similarity measures to quantify distances between probability density functions. In this regard, it aims to provide a framework for comparing
natural patterns in a statistical notation.  

This project is born out of my passion for statistics and I hope that it will be useful to
the people who share it with me.


## Tutorials 

 - [Introduction](https://github.com/HajkD/philentropy/blob/master/vignettes/Introduction.Rmd)
 - [Distances and Similarity Measures](https://github.com/HajkD/philentropy/blob/master/vignettes/Distances.Rmd)
 - [Information Theory](https://github.com/HajkD/philentropy/blob/master/vignettes/Information_Theory.Rmd)
 
 
### Installation

```r
# install.packages("devtools")
# install the current version of philentropy on your system
library(devtools)
install_github("HajkD/philentropy", build_vignettes = TRUE, dependencies = TRUE)
```

## Important Functions

### Distance Measures
* `distance()` : Implements 46 fundamental probability distance (or similarity) measures
* `getDistMethods()` : Get available method names for 'distance'
* `dist.diversity()` : Distance Diversity between Probability Density Functions
* `estimate.probability()` : Estimate Probability Vectors From Count Vectors

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



