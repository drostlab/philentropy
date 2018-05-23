---
title: 'Philentropy: Information Theory and Distance Quantification with R'
tags:
  - R
  - information theory
  - distance metrics
  - probability functions
  - divergence quantification
  - jensen-shannon divergence
authors:
  - name: Hajk-Georg Drost
    orcid: 0000-0002-1567-306X
    affiliation: "1"
affiliations:
 - name: The Sainsbury Laboratory, University of Cambridge, Bateman Street, Cambridge CB2 1LR, UK
   index: 1
date: 22 May 2018
bibliography: paper.bib
---

# Summary

Comparison is a fundamental method of scientific research leading to insights about the processes that generate similarity or dissimilarity. In statistical terms comparisons between probability functions are performed to infer connections, correlations, or relationships between objects or samples [@Cha2007]. Most quantification methods rely on distance or similarity measures, but the right choice for each individual application is not always clear and sometimes poorly explored.
The reason for this is partly that diverse measures are either implemented in different R packages with very different notations or are not implemented at all. Thus, a comprehensive framework implementing the most common similarity and distance measures using a uniform notation is still missing.
The R [@R2018] package ``Philentropy`` aims to fills this gap by implementing forty-six fundamental distance and similarity measures [@Cha2007] for comparing probability functions. These comparisons between probability functions have their foundations in a broad range of scientific disciplines from mathematics to ecology. The aim of this package is to provide a comprehensive and computationally optimized base framework for clustering, classification, statistical inference, goodness-of-fit, non-parametric statistics, information theory, and machine learning tasks that are based on comparing univariate or multivariate probability functions. All functions are written in C++ and are integrated into the R package using the Rcpp Application Programming Interface (API) [@Eddelbuettel2013].

Together, this framework allows building new similarity or distance based (statistical) models and algorithms in R which are computationally efficient and scalable. The comprehensive availability of diverse metrics and measures furthermore enables a systematic assessment of choosing the most optimal similarity or distance measure for individual applications in diverse scientific disciplines. 

The following probability distance/similarity and information theory measures are implemented in ``Philentropy``.

### Distance and Similarity Measures 

#### $L_p$ Minkowski Family
* Euclidean : $d = \sqrt{\sum_{i = 1}^N | P_i - Q_i |^2)}$
* Manhattan : $d = \sum_{i = 1}^N | P_i - Q_i |$
* Minkowski : $d = ( \sum_{i = 1}^N | P_i - Q_i |^p)^{1/p}$
* Chebyshev : $d = max | P_i - Q_i |$

#### $L_1$ Family
* Sorensen : $d = \frac{\sum_{i = 1}^N | P_i - Q_i |}{\sum_{i = 1}^N (P_i + Q_i)}$
* Gower : $d = \frac{1}{N} \dot \sum_{i = 1}^N | P_i - Q_i |$, where $N$ is the total number of elements $i$ in $P_i$ and $Q_i$
* Soergel : $d = \frac{\sum_{i = 1}^N | P_i - Q_i |}{\sum_{i = 1}^N max(P_i , Q_i)}$
* Kulczynski d : $d = \frac{\sum_{i = 1}^N | P_i - Q_i |}{\sum_{i = 1}^N min(P_i , Q_i)}$
* Canberra : $d = \frac{\sum_{i = 1}^N | P_i - Q_i |}{(P_i + Q_i)}$
* Lorentzian : $d = \sum_{i = 1}^N ln(1 + | P_i - Q_i |)$

#### Intersection Family
* Intersection : $s = \sum_{i = 1}^N min(P_i , Q_i)$
* Non-Intersection : $d = 1 - \sum_{i = 1}^N min(P_i , Q_i)$
* Wave Hedges : $d = \frac{\sum_{i = 1}^N | P_i - Q_i |}{max(P_i , Q_i)}$
* Czekanowski : $d = \frac{\sum_{i = 1}^N | P_i - Q_i |}{\sum_{i = 1}^N | P_i + Q_i |}$
* Motyka : $d = \frac{\sum_{i = 1}^N min(P_i , Q_i)}{(P_i + Q_i)}$
* Kulczynski s : $d = \frac{\sum_{i = 1}^N min(P_i , Q_i)}{\sum_{i = 1}^N | P_i - Q_i |}$
* Tanimoto : $d = \frac{\sum_{i = 1}^N (max(P_i , Q_i) - min(P_i , Q_i))}{\sum_{i = 1}^N max(P_i , Q_i)}$ ; equivalent to Soergel
* Ruzicka : $s = \frac{\sum_{i = 1}^N min(P_i , Q_i)}{\sum_{i = 1}^N max(P_i , Q_i)}$ ; equivalent to 1 - Tanimoto = 1 - Soergel 

#### Inner Product Family
* Inner Product : $s = \sum_{i = 1}^N P_i \dot Q_i$
* Harmonic mean : $s = 2 \cdot \frac{ \sum_{i = 1}^N P_i \cdot Q_i}{P_i + Q_i}$
* Cosine : $s = \frac{\sum_{i = 1}^N P_i \cdot Q_i}{\sqrt{\sum_{i = 1}^N P_i^2} \cdot \sqrt{\sum_{i = 1}^N Q_i^2}}$
* Kumar-Hassebrook (PCE) : $s = \frac{\sum_{i = 1}^N (P_i \cdot Q_i)}{(\sum_{i = 1}^N P_i^2 + \sum_{i = 1}^N Q_i^2 - \sum_{i = 1}^N (P_i \cdot Q_i))}$
* Jaccard : $d = 1 - \frac{\sum_{i = 1}^N P_i \cdot Q_i}{\sum_{i = 1}^N P_i^2 + \sum_{i = 1}^N Q_i^2 - \sum_{i = 1}^N P_i \cdot Q_i}$ ; equivalent to 1 - Kumar-Hassebrook
* Dice : $d = \frac{\sum_{i = 1}^N (P_i - Q_i)^2}{(\sum_{i = 1}^N P_i^2 + \sum_{i = 1}^N Q_i^2)}$

#### Squared-chord Family
* Fidelity : $s = \sum_{i = 1}^N \sqrt{P_i \cdot Q_i}$
* Bhattacharyya : $d = - ln \sum_{i = 1}^N \sqrt{P_i \cdot Q_i}$
* Hellinger : $d = 2 \cdot \sqrt{1 - \sum_{i = 1}^N \sqrt{P_i \cdot Q_i}}$
* Matusita : $d = \sqrt{2 - 2 \cdot \sum_{i = 1}^N \sqrt{P_i \cdot Q_i}}$
* Squared-chord : $d = \sum_{i = 1}^N ( \sqrt{P_i} - \sqrt{Q_i} )^2$

#### Squared $L_2$ family ($X^2$ squared family)
* Squared Euclidean : $d = \sum_{i = 1}^N ( P_i - Q_i )^2$
* Pearson $X^2$ : $d = \sum_{i = 1}^N ( \frac{(P_i - Q_i )^2}{Q_i} )$
* Neyman $X^2$ : $d = \sum_{i = 1}^N ( \frac{(P_i - Q_i )^2}{P_i} )$
* Squared $X^2$ : $d = \sum_{i = 1}^N ( \frac{(P_i - Q_i )^2}{(P_i + Q_i)} )$
* Probabilistic Symmetric $X^2$ : $d = 2 \cdot \sum_{i = 1}^N ( \frac{(P_i - Q_i )^2}{(P_i + Q_i)} )$
* Divergence : $X^2$ : $d = 2 \cdot \sum_{i = 1}^N ( \frac{(P_i - Q_i )^2}{(P_i + Q_i)^2} )$
* Clark : $d = \sqrt{\sum_{i = 1}^N (\frac{| P_i - Q_i |}{(P_i + Q_i)^2}}$
* Additive Symmetric $X^2$ : $d = \sum_{i = 1}^N ( \frac{((P_i - Q_i)^2 \cdot (P_i + Q_i))}{(P_i \cdot Q_i)} )$

#### Shannon's Entropy Family
* Kullback-Leibler : $d = \sum_{i = 1}^N P_i \cdot log(\frac{P_i}{Q_i})$
* Jeffreys : $d = \sum_{i = 1}^N (P_i - Q_i) \cdot log(\frac{P_i}{Q_i})$
* K divergence : $d = \sum_{i = 1}^N P_i \cdot log(\frac{2 \cdot P_i}{P_i + Q_i})$
* Topsoe : $d = \sum_{i = 1}^N ( P_i \cdot log(\frac{2 \cdot P_i}{P_i + Q_i}) ) + ( Q_i \cdot log(\frac{2 \cdot Q_i}{P_i + Q_i}) )$
* Jensen-Shannon :  $d = 0.5 \cdot ( \sum_{i = 1}^N P_i \cdot log(\frac{2 \cdot P_i}{P_i + Q_i}) + \sum_{i = 1}^N Q_i \cdot log(\frac{2 * Q_i}{P_i + Q_i}))$ 
* Jensen difference : $d = \sum_{i = 1}^N ( (\frac{P_i \cdot log(P_i) + Q_i \cdot log(Q_i)}{2}) - (\frac{P_i + Q_i}{2}) \cdot log(\frac{P_i + Q_i}{2}) )$

#### Combinations
* Taneja : $d = \sum_{i = 1}^N ( \frac{P_i + Q_i}{2}) \cdot log( \frac{P_i + Q_i}{( 2 \cdot \sqrt{P_i \cdot Q_i})} )$
* Kumar-Johnson : $d = \sum_{i = 1}^N \frac{(P_i^2 - Q_i^2)^2}{2 \cdot (P_i \cdot Q_i)^{\frac{3}{2}}}$
* Avg($L_1$, $L_n$) : $d = \frac{\sum_{i = 1}^N | P_i - Q_i| + max{ | P_i - Q_i |}}{2}$

__Note__: $d$ refers to distance measures, whereas $s$ denotes similarity measures.

### Information Theory Measures

- Shannon's Entropy H(X) : $H(X) = -\sum\limits_{i=1}^n P(x_i) \cdot log_b(P(x_i))$
- Shannon's Joint-Entropy H(X,Y) : $H(X,Y) = -\sum\limits_{i=1}^n\sum\limits_{j=1}^m P(x_i, y_j) \cdot log_b(P(x_i, y_j))$
- Shannon's Conditional-Entropy H(X | Y) : $H(Y|X) = \sum\limits_{i=1}^n\sum\limits_{j=1}^m P(x_i, y_j) \cdot log_b( \frac{P(x_i)}{P(x_i, y_j)})$
- Mutual Information I(X,Y) : $MI(X,Y) = \sum\limits_{i=1}^n\sum\limits_{j=1}^m P(x_i, y_j) \cdot log_b( \frac{P(x_i, y_j)}{( P(x_i) * P(y_j) )})$
- Kullback-Leibler Divergence : $KL(P || Q) = \sum\limits_{i=1}^n P(p_i) \cdot log_2(\frac{P(p_i) }{P(q_i)}) = H(P, Q) - H(P)$
- Jensen-Shannon Divergence : $JSD(P || Q) = 0.5 * (KL(P || R) + KL(Q || R))$
- Generalized Jensen-Shannon Divergence : $gJSD_{\pi_1,...,\pi_n}(P_1, ..., P_n) = H(\sum_{i = 1}^n \pi_i \cdot P_i) - \sum_{i = 1}^n \pi_i \cdot H(P_i)$

### Example 

The following example shall illustrate how to use the ``Philentropy`` package.

```{R}
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

#### Example: Jensen-Shannon Divergence

This function computes the `Jensen-Shannon Divergence` `JSD(P || Q)` between two probability distributions `P` and `Q` with equal weights `π_1` = `π_2` = 1/2.

The Jensen-Shannon Divergence `JSD(P || Q)` between two probability distributions P and Q is defined as:

$$JSD(P || Q) = 0.5 * (KL(P || R) + KL(Q || R))$$

where `R = 0.5 * (P + Q)` denotes the mid-point of the probability vectors `P` and `Q`, and `KL(P || R)`, `KL(Q || R)` denote the `Kullback-Leibler Divergence` of `P` and `R`, as well as `Q` and `R`.

```r
# Jensen-Shannon Divergence between P and Q
P <- 1:10/sum(1:10)
Q <- 20:29/sum(20:29)
x <- rbind(P,Q)

# Jensen-Shannon Divergence between P and Q using different log bases
distance(x, method = "jensen-shannon", unit = "log2") # Default
distance(x, method = "jensen-shannon", unit = "log")
distance(x, method = "jensen-shannon", unit = "log10")
```

```
# Metric: 'jensen-shannon' using unit: 'log2'.
jensen-shannon 
    0.03792749 
# Metric: 'jensen-shannon' using unit: 'log'.
jensen-shannon 
    0.02628933 
# Metric: 'jensen-shannon' using unit: 'log10'.
jensen-shannon 
    0.01141731 
```

Alternatively, users can specify count data.

```{R}
# Jensen-Shannon Divergence Divergence between count vectors P_count and Q_count
P_count <- 1:10
Q_count <- 20:29
x_count <- rbind(P_count, Q_count)

distance(x_count, method = "jensen-shannon", est.prob = "empirical")
```

```
# Metric: 'jensen-shannon' using unit: 'log'.
jensen-shannon 
    0.02628933
```

Or users can compute distances based on a probability matrix

```r
# Example: Distance Matrix using JSD-Distance
Prob <- rbind(1:10/sum(1:10), 20:29/sum(20:29), 30:39/sum(30:39))

# compute the JSD matrix of a given probability matrix
JSDMatrix <- distance(Prob, method = "jensen-shannon", unit = "log2")

JSDMatrix
```

```
# Metric: 'jensen-shannon' using unit: 'log2'.
           v1           v2           v3
v1 0.00000000 0.0379274917 0.0435852218
v2 0.03792749 0.0000000000 0.0002120578
v3 0.04358522 0.0002120578 0.0000000000
```

#### Properties of the `Jensen-Shannon Divergence`:

- JSD is non-negative.

- JSD is a symmetric measure JSD(P || Q) = JSD(Q || P).

- JSD = 0, if and only if P = Q.

``Philentropy`` already enabled the robust comparison of similarity measures in
analogy-based software effort estimation [@Phannachitta2017] as well as in evolutionary transcriptomics applications [@Drost2018]. The package aims to assist efforts 
to determine the optimal similarity or distance measures when developing new (statistical) models or algorithms. In addition, ``Philentropy`` is efficiently implemented to be applicable to large-scale datasets that were previously inaccessible using other R packages. The software is open source and currently available on GitHub (https://github.com/HajkD/philentropy) and CRAN (https://cran.r-project.org/web/packages/philentropy/index.html).


# Acknowledgements

I would like to thank Jerzy Paszkowski for providing me an inspiring scientific environment
and for supporting my projects.

# References
