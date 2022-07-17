philentropy
===========

[![Travis-CI Build Status](https://travis-ci.org/HajkD/philentropy.svg?branch=master)](https://travis-ci.org/HajkD/philentropy)  [![status](http://joss.theoj.org/papers/cad5ffc246ce197b06ccad1af7d2932a/status.svg)](https://joss.theoj.org/papers/10.21105/joss.00765)
[![rstudio mirror downloads](http://cranlogs.r-pkg.org/badges/philentropy)](https://github.com/r-hub/cranlogs.app)
[![rstudio mirror downloads](http://cranlogs.r-pkg.org/badges/grand-total/philentropy)](https://github.com/r-hub/cranlogs.app)


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
# install philentropy version 0.6.0 from CRAN
install.packages("philentropy")
```

### Citation

__I am developing `philentropy` in my spare time and would be very grateful if you would consider citing the following paper in case `philentropy` was useful for your own research. I plan on maintaining and extending the `philentropy` functionality and usability in the next years and require citations to back up these efforts. Many thanks in advance :)__

> HG Drost, (2018). __Philentropy: Information Theory and Distance Quantification with R__. _Journal of Open Source Software_, 3(26), 765. https://doi.org/10.21105/joss.00765

## Tutorials 

 - [Introduction to the philentropy package](https://drostlab.github.io/philentropy/articles/Introduction.html)
 - [Distance and Similarity Measures implemented in philentropy](https://drostlab.github.io/philentropy/articles/Distances.html)
 - [Information Theory Metrics implemented in philentropy](https://drostlab.github.io/philentropy/articles/Information_Theory.html)
- [Comparing many probability density functions](https://drostlab.github.io/philentropy/articles/Many_Distances.html)

## Examples

```r
library(philentropy)
# retrieve available distance metrics
philentropy::getDistMethods()
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
philentropy::distance(x, method = "jensen-shannon")
```

```
jensen-shannon using unit 'log'.
jensen-shannon 
    0.02628933
```

Alternatively, users can also retrieve values from all available distance/similarity metrics
using `philentropy::dist.diversity()`:


```r
philentropy::dist.diversity(x, p = 2, unit = "log2")
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

The current status of the package as well as a detailed history of the functionality of each version of `philentropy` can be found in the [NEWS](https://drostlab.github.io/philentropy/news/index.html) section.

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

> - __An atlas of gene regulatory elements in adult mouse cerebrum__ YE Li, S Preissl, X Hou, Z Zhang, K Zhang et al.- __Nature__, 2021
>
> - __Convergent somatic mutations in metabolism genes in chronic liver disease__ S Ng, F Rouhani, S Brunner, N Brzozowska et al. __Nature__, 2021
>
> - __Antigen dominance hierarchies shape TCF1+ progenitor CD8 T cell phenotypes in tumors__ ML Burger, AM Cruz, GE Crossland et al. - __Cell__, 2021
>
> - __High-content single-cell combinatorial indexing__ R Mulqueen et al. - __Nature Biotechnology__, 2021 
>
> - __Extinction at the end-Cretaceous and the origin of modern Neotropical rainforests__ MR Carvalho, C Jaramillo et al. - __Science__, 2021
>
> - __HERMES: a molecular-formula-oriented method to target the metabolome__
R Giné, J Capellades, JM Badia et al. - __Nature Methods__, 2021
>
> - __The genetic architecture of temperature adaptation is shaped by population ancestry and not by selection regime__ KA Otte, V Nolte, F Mallard et al. - __Genome Biology__, 2021
>
> - __The Tug1 lncRNA locus is essential for male fertility__ JP Lewandowski et al. - __Genome Biology__, 2020
> 
> - __Resolving the structure of phage–bacteria interactions in the context of natural diversity__ KM Kauffman, WK Chang, JM Brown et al. - __Nature Communications__, 2022
>
> - __Gut microbiome-mediated metabolism effects on immunity in rural and urban African populations__
M Stražar, GS Temba, H Vlamakis et al. - __Nature Communications__, 2021
>
> - __Aging, inflammation and DNA damage in the somatic testicular niche with idiopathic germ cell aplasia__ M Alfano, AS Tascini, F Pederzoli et al. - __Nature communications__, 2021
>
> - __Single cell census of human kidney organoids shows reproducibility and diminished off-target cells after transplantation__ A Subramanian et al. - __Nature Communications__, 2019
>
> - __Different languages, similar encoding efficiency: Comparable information rates across the human communicative niche__
C Coupé, YM Oh, D Dediu, F Pellegrino - __Science Advances__, 2019
>
> - __Single-cell deletion analyses show control of pro–T cell developmental speed and pathways by Tcf7, Spi1, Gata3, Bcl11a, Erg, and Bcl11b__ W Zhou, F Gao, M Romero-Wolf, S Jo, EV Rothenberg - __Science Immunology__, 2022
>
> - __Large-scale chromatin reorganization reactivates placenta-specific genes that drive cellular aging__ Z Liu, Q Ji, J Ren, P Yan, Z Wu, S Wang, L Sun, Z Wang et al - __Developmental Cell__, 2022
>
> - __Direct epitranscriptomic regulation of mammalian translation initiation through N4-acetylcytidine__ D Arango, D Sturgill, R Yang, T Kanai, P Bauer et al. - __Molecular Cell__, 2022
>
> - __Loss of adaptive capacity in asthmatic patients revealed by biomarker fluctuation dynamics after rhinovirus challenge__ A Sinha et al. - __eLife__, 2019
>
> - __Sex and hatching order modulate the association between MHC‐II diversity and fitness in early‐life stages of a wild seabird__
M Pineaux et al - __Molecular Ecology__, 2020
>
> - __How the Choice of Distance Measure Influences the Detection of Prior-Data Conflict__
K Lek, R Van De Schoot - __Entropy__, 2019
>
> - __Differential variation analysis enables detection of tumor heterogeneity using single-cell RNA-sequencing data__
EF Davis-Marcisak, TD Sherman et al. - __Cancer research__, 2019
>
> - __Multi-Omics Investigation of Innate Navitoclax Resistance in Triple-Negative Breast Cancer Cells__ M Marczyk et al. - __Cancers__, 2020
>
> - __Impact of Gut Microbiome on Hypertensive Patients with Low-Salt Intake: Shika Study Results__ S Nagase et al. - __Frontiers in Medicine__, 2020
>
> - __Combined TCR Repertoire Profiles and Blood Cell Phenotypes Predict Melanoma Patient Response to Personalized Neoantigen Therapy plus Anti-PD-1__ A Poran et al. - __Cell Reports Medicine__, 2020
>
> - __Identification of a glioma functional network from gene fitness data using machine learning__ C Xiang, X Liu, D Zhou, Y Zhou, X Wang, F Chen - __Journal of Cellular and Molecular Medicine__, 2022
>
> - __Prediction of New Risk Genes and Potential Drugs for Rheumatoid Arthritis from Multiomics Data__ AM Birga, L Ren, H Luo, Y Zhang, J Huang - __Computational and Mathematical Methods in Medicine__, 2022
>
> - __Phenotyping of acute and persistent COVID-19 features in the outpatient setting: exploratory analysis of an international cross-sectional online survey__ S Sahanic, P Tymoszuk, D Ausserhofer et al. - __medRxiv__, 2021
>
> - __A two-part evaluation approach for measuring the usability and user experience of an Augmented Reality-based assistance system to support the temporal coordination of spatially dispersed teams__ L Thomaschewski, B Weyers, A Kluge - __Cognitive Systems Research__, 2021
>
> - __SEDE-GPS: socio-economic data enrichment based on GPS information__
T Sperlea, S Füser, J Boenigk, D Heider - __BMC bioinformatics__, 2018
>
> - __Spatial and molecular anatomy of germ layers in the gastrulating primate embryo__ G Cui, S Feng, Y Yan, L Wang, X He, X Li, et al. - __bioRxiv__, 2022
>
> - __Evacuees and Migrants Exhibit Different Migration Systems after the Great East Japan Earthquake and Tsunami__
M Hauer, S Holloway, T Oda – 2019
>
> - __Robust comparison of similarity measures in analogy based software effort estimation__
P Phannachitta - __11th International Conference on Software__, 2017
>
> - __RUNIMC - An R-based package for imaging mass cytometry data analysis and pipeline validation__ L Dolcetti, PR Barber, G Weitsman, S Thavarajet al. - __bioRxiv__, 2021
>
> - __Expression variation analysis for tumor heterogeneity in single-cell RNA-sequencing data__
EF Davis-Marcisak, P Orugunta et al. - __BioRxiv__, 2018
>
> - __Concept acquisition and improved in-database similarity analysis for medical data__
I Wiese, N Sarna, L Wiese, A Tashkandi, U Sax - __Distributed and Parallel Databases__, 2019
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
>
> - __Multi-classifier majority voting analyses in provenance studies on iron artefacts__
G Żabiński et al. - __Journal of Archaeological Science__, 2020
>
> - __Identifying inhibitors of epithelial–mesenchymal plasticity using a network topology-based approach__
K Hari et al. - __NPJ systems biology and applications__, 2020
>
> - __Genetic differentiation and intrinsic genomic features explain variation in recombination hotspots among cocoa tree populations__
EJ Schwarzkopf et al. - __BMC Genomics__, 2020
>
> - __Enhancing Card Sorting Dendrograms through the Holistic Analysis of Distance Methods and Linkage Criteria.__ JA Macías - __Journal of Usability Studies__, 2021
>
> - __Pattern-based identification and mapping of landscape types using multi-thematic data__ J Nowosad, TF Stepinski - __International Journal of Geographical Information__, 2021
>
> - __Motif Analysis in k-mer Networks: An Approach towards Understanding SARS-CoV-2 Geographical Shifts__
S Biswas, S Saha, S Bandyopadhyay, M Bhattacharyya - __bioRxiv__, 2020
>
> - __Motif: an open-source R tool for pattern-based spatial analysis__ J Nowosad - __Landscape Ecology__, 2021
>
> - __New effective spectral matching measures for hyperspectral data analysis__ C Kumar, S Chatterjee, T Oommen, A Guha - __International Journal of Remote Sensing__, 2021
>
> - __Innovative activity of Polish enterprises–a strategic aspect. The similarity of NACE divisions__ E Bielińska-Dusza, M Hamerska - __Journal of Entrepreneurship, Management and innovation__, 2021
>
> - __Multi-classifier majority voting analyses in provenance studies on iron artefacts__ G Żabiński, J Gramacki et al.- __Journal of Archaeological Science__, 2020
>
> - __Unraveling the record of a tropical continental Cretaceous-Paleogene boundary in northern Colombia, South America__ F de la Parra, C Jaramillo, P Kaskes et al. - __Journal of South American Earth Sciences__, 2022
>
> - __A roadmap to reconstructing muscle architecture from CT data__
J Katzke, P Puchenkov, H Stark, EP Economo - __Integrative Organismal Biology__, 2022
> 
> - __Pandemonium: a clustering tool to partition parameter space—application to the B anomalies__ U Laa, G Valencia - __The European Physical Journal Plus__, 2022
>
> - __Identification of a glioma functional network from gene fitness data using machine learning__ C Xiang, X Liu, D Zhou, Y Zhou, X Wang, F Chen - __Journal of Cellular and Molecular Medicine__, 2022
>
> - __Cross compatibility in intraspecific and interspecific hybridization in yam (Dioscorea spp.)__ JM Mondo, PA Agre, A Edemodu et al. - __Scientific reports__, 2022
>
> - __A Modular and Expandable Ecosystem for Metabolomics Data Annotation in R__ J Rainer, A Vicini, L Salzer, J Stanstrup et al. - __Metabolites__, 2022
>
> - __Single-Cell Transcriptome Integration Analysis Reveals the Correlation Between Mesenchymal Stromal Cells and Fibroblasts__ C Fan, M Liao, L Xie, L Huang, S Lv, S Cai et al. - __Frontiers in genetics__, 2022
>
> - __Phenotypic regionalization of the vertebral column in the thorny skate Amblyraja radiata: Stability and variation__ F Berio, Y Bayle, C Riley, O Larouche, R Cloutier - __Journal of Anatomy__, 2022
>
> - __Community assembly during vegetation succession after metal mining is driven by multiple processes with temporal variation__ T Li, H Yang, X Yang, Z Guo, D Fu, C Liu, S Li et al. - __Ecology and evolution__, 2022
>
> - __Integrative Organismal Biology__ J Katzke, P Puchenkov, H Stark, __EP Economo__ - 2022
>
> - __Optimizing use of US Ex-PVP inbred lines for enhancing agronomic performance of tropical Striga resistant maize inbred lines__ ARS Maazou, M Gedil, VO Adetimirin et al. - __BMC Plant Biology__, 2022



## Discussions and Bug Reports

I would be very happy to learn more about potential improvements of the concepts and functions
provided in this package.

Furthermore, in case you find some bugs or need additional (more flexible) functionality of parts
of this package, please let me know:

https://github.com/drostlab/philentropy/issues

or find me on [twitter: HajkDrost](https://twitter.com/hajkdrost) 



