philentropy
===========

[![Travis-CI Build Status](https://travis-ci.org/HajkD/philentropy.svg?branch=master)](https://travis-ci.org/HajkD/philentropy) [![rpackages.io rank](https://www.rpackages.io/badge/philentropy.svg)](https://www.rpackages.io/package/philentropy)


### Similarity and Distance Quantification between Probability Functions

> Describe and understand the world through data.

Data collection and data comparison are the foundations of scientific research.
_Mathematics_ provides the abstract framework to describe patterns we observe in nature and _Statistics_ provides the
framework to quantify the uncertainty of these patterns. In statistics, natural patterns
are described in form of probability distributions which either follow a fixed pattern (parametric distributions) or more dynamic patterns (non-parametric distributions).

The `philentropy` package implements fundamental distance and similarity measures to quantify distances between probability density functions as well as traditional information theory measures. In this regard, it aims to provide a framework for comparing
natural patterns in a statistical notation.  

This project is born out of my passion for statistics and I hope that it will be useful to
the people who share it with me.

### Installation
```r
# install philentropy version 0.1.0 from CRAN
install.packages("philentropy")
```

## Tutorials 

 - [Introduction to the philentropy package](https://hajkd.github.io/philentropy/articles/Introduction.html)
 - [Distance and Similarity Measures implemented in philentropy](https://hajkd.github.io/philentropy/articles/Distances.html)
 - [Information Theory Metrics implemented in philentropy](https://hajkd.github.io/philentropy/articles/Information_Theory.html)

## Examples

```r
library(philentropy)
# retrieve available distance metrics
getDistMethods()
```

```
 [1] "euclidean"         "manhattan"         "minkowski"        
 [4] "chebyshev"         "sorensen"          "gower"            
 [7] "soergel"           "kulczynski_d"      "canberra"         
[10] "lorentzian"        "intersection"      "non-intersection" 
[13] "wavehedges"        "czekanowski"       "motyka"           
[16] "kulczynski_s"      "tanimoto"          "ruzicka"          
[19] "inner_product"     "harmonic_mean"     "cosine"           
[22] "hassebrook"        "jaccard"           "dice"             
[25] "fidelity"          "bhattacharyya"     "hellinger"        
[28] "matusita"          "squared_chord"     "squared_euclidean"
[31] "pearson"           "neyman"            "squared_chi"      
[34] "prob_symm"         "divergence"        "clark"            
[37] "additive_symm"     "kullback-leibler"  "jeffreys"         
[40] "k_divergence"      "topsoe"            "jensen-shannon"   
[43] "jensen_difference" "taneja"            "kumar-johnson"    
[46] "avg"
```

```r
# define a probability density function P
P <- 1:10/sum(1:10)
# define a probability density function Q
Q <- 20:29/sum(20:29)

# combine P and Q as matrix object
x <- rbind(P,Q)

# compute the jensen-shannon distance between
# probability density functions P and Q
distance(x, method = "jensen-shannon")
```

```
jensen-shannon using unit 'log'.
jensen-shannon 
    0.02628933
```
 
### Install Developer Version
```r
# install.packages("devtools")
# install the current version of philentropy on your system
library(devtools)
install_github("HajkD/philentropy", build_vignettes = TRUE, dependencies = TRUE)
```

### NEWS

The current status of the package as well as a detailed history of the functionality of each version of `philentropy` can be found in the [NEWS](https://hajkd.github.io/philentropy/news/index.html) section.

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


## Discussions and Bug Reports

I would be very happy to learn more about potential improvements of the concepts and functions
provided in this package.

Furthermore, in case you find some bugs or need additional (more flexible) functionality of parts
of this package, please let me know:

https://github.com/HajkD/philentropy/issues

or find me on [twitter: HajkDrost](https://twitter.com/hajkdrost) 



