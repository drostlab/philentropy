philentropy
===========

[![Travis-CI Build Status](https://travis-ci.org/HajkD/philentropy.svg?branch=master)](https://travis-ci.org/HajkD/philentropy)  [![status](http://joss.theoj.org/papers/cad5ffc246ce197b06ccad1af7d2932a/status.svg)](http://joss.theoj.org/papers/cad5ffc246ce197b06ccad1af7d2932a)
[![rstudio mirror downloads](http://cranlogs.r-pkg.org/badges/philentropy)](https://github.com/metacran/cranlogs.app)
[![rstudio mirror downloads](http://cranlogs.r-pkg.org/badges/grand-total/philentropy)](https://github.com/metacran/cranlogs.app)


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
# install philentropy version 0.3.0 from CRAN
install.packages("philentropy")
```

### Citation

__I am developing `philentropy` in my spare time and would be very grateful if you would consider citing the following paper in case `philentropy` was useful for your own research. I plan on maintaining and extending the `philentropy` functionality and usability in the next years and require citations to back up these efforts. Many thanks in advance :)__

> HG Drost, (2018). __Philentropy: Information Theory and Distance Quantification with R__. _Journal of Open Source Software_, 3(26), 765. https://doi.org/10.21105/joss.00765

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

Alternatively, users can also retrieve values from all available distance/similarity metrics
using `dist.diversity()`:


```r
dist.diversity(x, p = 2, unit = "log2")
```

```
        euclidean         manhattan 
       0.12807130        0.35250464 
        minkowski         chebyshev 
       0.12807130        0.06345083 
         sorensen             gower 
       0.17625232        0.03525046 
          soergel      kulczynski_d 
       0.29968454        0.42792793 
         canberra        lorentzian 
       2.09927095        0.49712136 
     intersection  non-intersection 
       0.82374768        0.17625232 
       wavehedges       czekanowski 
       3.16657887        0.17625232 
           motyka      kulczynski_s 
       0.58812616        2.33684211 
         tanimoto           ruzicka 
       0.29968454        0.70031546 
    inner_product     harmonic_mean 
       0.10612245        0.94948528 
           cosine        hassebrook 
       0.93427641        0.86613103 
          jaccard              dice 
       0.13386897        0.07173611 
         fidelity     bhattacharyya 
       0.97312397        0.03930448 
        hellinger          matusita 
       0.32787819        0.23184489 
    squared_chord squared_euclidean 
       0.05375205        0.01640226 
          pearson            neyman 
       0.16814418        0.36742465 
      squared_chi         prob_symm 
       0.10102943        0.20205886 
       divergence             clark 
       1.49843905        0.86557468 
    additive_symm  kullback-leibler 
       0.53556883        0.13926288 
         jeffreys      k_divergence 
       0.31761069        0.04216273 
           topsoe    jensen-shannon 
       0.07585498        0.03792749 
jensen_difference            taneja 
       0.03792749        0.04147518 
    kumar-johnson               avg 
       0.62779644        0.20797774
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

### Studies that successfully applied the `philentropy` package

> - __Single cell census of human kidney organoids shows reproducibility and diminished off-target cells after transplantation__ A Subramanian et al. - __Nature Communications__, 2019
>
> - __Different languages, similar encoding efficiency: Comparable information rates across the human communicative niche__
C Coupé, YM Oh, D Dediu, F Pellegrino - __Science Advances__, 2019
>
> - __Loss of adaptive capacity in asthmatic patients revealed by biomarker fluctuation dynamics after rhinovirus challenge__ A Sinha et al. - __eLife__, 2019
>
> - __Evacuees and Migrants Exhibit Different Migration Systems after the Great East Japan Earthquake and Tsunami__
M Hauer, S Holloway, T Oda – 2019
>
> - __Robust comparison of similarity measures in analogy based software effort estimation__
P Phannachitta - __11th International Conference on Software__, 2017
>
> - __Expression variation analysis for tumor heterogeneity in single-cell RNA-sequencing data__
EF Davis-Marcisak, P Orugunta et al. - __BioRxiv__, 2018
>
> - __SEDE-GPS: socio-economic data enrichment based on GPS information__
T Sperlea, S Füser, J Boenigk, D Heider - __BMC bioinformatics__, 2018
>
> - __How the Choice of Distance Measure Influences the Detection of Prior-Data Conflict__
K Lek, R Van De Schoot - __Entropy__, 2019
>
> - __Concept acquisition and improved in-database similarity analysis for medical data__
I Wiese, N Sarna, L Wiese, A Tashkandi, U Sax - __Distributed and Parallel Databases__, 2019
>
> - __Differential variation analysis enables detection of tumor heterogeneity using single-cell RNA-sequencing data__
EF Davis-Marcisak, TD Sherman et al. - __Cancer research__, 2019
>
> - __Dynamics of Vaginal and Rectal Microbiota over Several Menstrual Cycles in Female Cynomolgus Macaques__
MT Nugeyre, N Tchitchek, C Adapen et al. - __Frontiers in Cellular and Infection Microbiology__, 2019
>
> - __Inferring the quasipotential landscape of microbial ecosystems with topological data analysis__
WK Chang, L Kelly - __BioRxiv__, 2019
>
> - __Shifts in the nasal microbiota of swine in response to different dosing regimens of oxytetracycline administration__
KT Mou, HK Allen, DP Alt, J Trachsel et al. - __Veterinary microbiology__, 2019
>
> - __The Patchy Distribution of Restriction–Modification System Genes and the Conservation of Orphan Methyltransferases in Halobacteria__
MS Fullmer, M Ouellette, AS Louyakis et al. - __Genes__, 2019
>
> - __Genetic differentiation and intrinsic genomic features explain variation in recombination hotspots among cocoa tree populations__
EJ Schwarzkopf, JC Motamayor, OE Cornejo - __BioRxiv__, 2019
>
> - __Metastable regimes and tipping points of biochemical networks with potential applications in precision medicine__
SS Samal, J Krishnan, AH Esfahani et al. - __Reasoning for Systems Biology and Medicine__, 2019
>
> - __Genome‐wide characterization and developmental expression profiling of long non‐coding RNAs in Sogatella furcifera__
ZX Chang, OE Ajayi, DY Guo, QF Wu - __Insect science__, 2019
>
> - __Loss of adaptive capacity in asthmatics revealed by biomarker fluctuation dynamics upon experimental rhinovirus challenge__
A Sinha, R Lutter, B Xu, T Dekker, B Dierdorp et al. - __BioRxiv__, 2019
>
> - __Development of a simulation system for modeling the stock market to study its characteristics__
P Mariya – 2018
>
> - __The Tug1 Locus is Essential for Male Fertility__
JP Lewandowski, G Dumbović, AR Watson, T Hwang et al. - __BioRxiv__, 2019
>
> - __Microbiotyping the sinonasal microbiome__
A Bassiouni, S Paramasivan, A Shiffer et al. - __BioRxiv__, 2019
>
> - __Critical search: A procedure for guided reading in large-scale textual corpora__
J Guldi - __Journal of Cultural Analytics__, 2018
>
> - __A Bibliography of Publications about the R, S, and S-Plus Statistics Programming Languages__
NHF Beebe – 2019
>
> - __Improved state change estimation in dynamic functional connectivity using hidden semi-Markov models__
H Shappell, BS Caffo, JJ Pekar, MA Lindquist - __NeuroImage__, 2019
>
> - __A Smart Recommender Based on Hybrid Learning Methods for Personal Well-Being Services__
RM Nouh, HH Lee, WJ Lee, JD Lee - __Sensors__, 2019
>
> - __Cognitive Structural Accuracy__
V Frenz – 2019
>
> - __Kidney organoid reproducibility across multiple human iPSC lines and diminished off target cells after transplantation revealed by single cell transcriptomics__
A Subramanian, EH Sidhom, M Emani et al. - __BioRxiv__, 2019
> - __Multi-classifier majority voting analyses in provenance studies on iron artefacts__
G Żabiński et al. - __Journal of Archaeological Science__, 2020






## Discussions and Bug Reports

I would be very happy to learn more about potential improvements of the concepts and functions
provided in this package.

Furthermore, in case you find some bugs or need additional (more flexible) functionality of parts
of this package, please let me know:

https://github.com/HajkD/philentropy/issues

or find me on [twitter: HajkDrost](https://twitter.com/hajkdrost) 



