# philentropy <sub><sup>— Information Theory and Distance Quantification with R</sup></sub>  

[![CRAN](https://img.shields.io/cran/v/philentropy)](https://cran.r-project.org/package=philentropy)
 [![status](http://joss.theoj.org/papers/cad5ffc246ce197b06ccad1af7d2932a/status.svg)](https://joss.theoj.org/papers/10.21105/joss.00765) 
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/philentropy)](https://github.com/r-hub/cranlogs.app)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/grand-total/philentropy)](https://github.com/r-hub/cranlogs.app)


## 🧭 Similarity and Distance Quantification between Probability Functions

> _Describe and understand the world through data._

Data collection and data comparison are the foundations of scientific research.  
_Mathematics_ provides the abstract framework to describe patterns we observe in nature and _Statistics_ provides the
framework to quantify the uncertainty of these patterns.  

In statistics, natural patterns are described in the form of probability distributions that either follow fixed patterns (parametric distributions) or more dynamic ones (non-parametric distributions).

The `philentropy` package implements fundamental distance and similarity measures to quantify distances between probability density functions as well as traditional information theory measures.  
In this regard, it aims to provide a framework for comparing natural patterns in a statistical notation.

> 🧡 This project is born out of my passion for statistics and I hope it will be useful to those who share it with me.

---

## ⚙️ Installation

```r
# install philentropy version 0.10.0 from CRAN
install.packages("philentropy")
```

Or get the latest developer version:

```r
# install.packages("devtools")
library(devtools)
install_github("HajkD/philentropy", build_vignettes = TRUE, dependencies = TRUE)
```

---

## 🧾 Citation

> HG Drost (2018).  
> __Philentropy: Information Theory and Distance Quantification with R__.  
> _Journal of Open Source Software_, 3(26), 765.  
> https://doi.org/10.21105/joss.00765

> 🪶 *I am developing `philentropy` in my spare time and would be very grateful if you would consider citing the paper above if it was useful for your research. These citations help me continue maintaining and extending the package.*

---

## 🧩 Quick Start

```r
library(philentropy)

P <- c(0.1, 0.2, 0.7)
Q <- c(0.2, 0.2, 0.6)

distance(rbind(P, Q), method = "jensen-shannon")
```

```
jensen-shannon using unit 'log'.
jensen-shannon 
    0.02628933
```

> 💡 **Tip:** Got a large matrix (rows = samples, cols = features)?  
> Use `distance(X, method="cosine", mute.message=TRUE)` to compute the full pairwise matrix quickly and quietly.

---

## 📘 Tutorials

- [Introduction to the philentropy package](https://drostlab.github.io/philentropy/articles/Introduction.html)  
- [Distance and Similarity Measures implemented in philentropy](https://drostlab.github.io/philentropy/articles/Distances.html)  
- [Information Theory Metrics implemented in philentropy](https://drostlab.github.io/philentropy/articles/Information_Theory.html)  
- [Comparing many probability density functions](https://drostlab.github.io/philentropy/articles/Many_Distances.html)

---

## 🧪 When should I use which distance?

| Goal | Recommended Methods |
|------|---------------------|
| 🔁 Clustering / similarity | `cosine`, `correlation`, `euclidean` |
| 📊 Probability or compositional data | `jensen-shannon`, `hellinger`, `kullback-leibler` |
| 🧬 Sparse counts / binary | `canberra`, `jaccard`, `sorensen` |
| ⚖️ Scale-invariant | `manhattan`, `chebyshev` |

> Run `getDistMethods()` to explore all 45+ implemented measures.

---

## 🧮 Examples

```r
library(philentropy)
philentropy::getDistMethods()
```

```
[1] "euclidean"         "manhattan"         "minkowski"         "chebyshev"         "sorensen"         
[6] "gower"             "soergel"           "kulczynski_d"      "canberra"          "lorentzian"       
[11] "intersection"      "non-intersection"  "wavehedges"        "czekanowski"       "motyka"           
[16] "kulczynski_s"      "tanimoto"          "ruzicka"           "inner_product"     "harmonic_mean"    
[21] "cosine"            "hassebrook"        "jaccard"           "dice"              "fidelity"         
[26] "bhattacharyya"     "hellinger"         "matusita"          "squared_chord"     "squared_euclidean"
[31] "pearson"           "neyman"            "squared_chi"       "prob_symm"         "divergence"       
[36] "clark"             "additive_symm"     "kullback-leibler"  "jeffreys"          "k_divergence"     
[41] "topsoe"            "jensen-shannon"    "jensen_difference" "taneja"            "kumar-johnson"    
[46] "avg"
```

```r
# define probability density functions P and Q
P <- 1:10/sum(1:10)
Q <- 20:29/sum(20:29)

x <- rbind(P, Q)
philentropy::distance(x, method = "jensen-shannon")
```

```
jensen-shannon using unit 'log'.
jensen-shannon 
    0.02628933
```

Alternatively, compute all available distances:

```r
philentropy::dist.diversity(x, p = 2, unit = "log2")
```

---

## 🌟 Papers using philentropy (highlights)

> _Flagship examples with top venues. Click to expand full lists._

<details>
<summary><b>Nature / Cell / Science</b></summary>

- __A transcriptomic hourglass in brown algae__  
  JS Lotharukpong, M Zheng, R Luthringer et al. – __Nature__, 2024  
- __Annelid functional genomics reveal the origins of bilaterian life cycles__  
  FM Martín-Zamora, Y Liang, K Guynes et al. – __Nature__, 2023  
- __An atlas of gene regulatory elements in adult mouse cerebrum__  
  YE Li, S Preissl, X Hou, Z Zhang, K Zhang et al. – __Nature__, 2021  
- __Convergent somatic mutations in metabolism genes in chronic liver disease__  
  S Ng, F Rouhani, S Brunner, N Brzozowska et al. – __Nature__, 2021  
- __Antigen dominance hierarchies shape TCF1+ progenitor CD8 T cell phenotypes in tumors__  
  ML Burger, AM Cruz, GE Crossland et al. – __Cell__, 2021  
- __A comparative atlas of single-cell chromatin accessibility in the human brain__  
  YE Li, S Preissl, M Miller, ND Johnson, Z Wang et al. – __Science__, 2023
</details>

<details>
<summary><b>Nature Methods / Nat Comms / Cell family</b></summary>

- __sciCSR infers B cell state transition and predicts class-switch recombination dynamics using scRNA-seq__  
  JCF Ng, G Montamat Garcia, AT Stewart et al. – __Nature Methods__, 2024  
- __Decoding the gene regulatory network of endosperm differentiation in maize__  
  Y Yuan, Q Huo, Z Zhang, Q Wang, J Wang et al. – __Nature Communications__, 2024  
- __Population structure in a fungal human pathogen is potentially linked to pathogenicity__  
  EA Hatmaker, AE Barber, MT Drott et al. – __Nature Communications__, 2025  
- __Pan-cancer human brain metastases atlas at single-cell resolution__  
  X Xing, J Zhong, J Biermann, H Duan, X Zhang et al. – __Cancer Cell__, 2025  
- __Gene module reconstruction identifies cellular differentiation processes and the regulatory logic of specialized secretion in zebrafish__  
  Y Wang, J Liu, LY Du, JL Wyss, JA Farrell, AF Schier – __Developmental Cell__, 2025  
</details>

<details>
<summary><b>Other disciplines (selected)</b></summary>

- __Staphylococci in high resolution: Capturing diversity within the human nasal microbiota__  
  AC Ingham, DYK Ng, S Iversen, CM Liu et al. – __Cell Reports__, 2025  
- __The power of visualizing distributional differences: formal graphical n-sample tests__  
  K Konstantinou, T Mrkvička, M Myllymäki – __Computational Statistics__, 2025  
- __Plant species as ecological engineers of microtopography in a temperate sedge-grass marsh__  
  J Dušek, J Novotný, B Navrátilová et al. – __Scientific Reports__, 2025  
- __Resolution of MALDI-TOF vs WGS for Bacillus identification (NASA JSC)__  
  F Mazhari, AB Regberg, CL Castro, MG LaMontagne – __Frontiers in Microbiology__, 2025  
- __Every Hue Has Its Fan Club: Diverse Patterns of Color-Dependent Flower Visitation across Lepidoptera__  
  D Kutcherov, EL Westerman – __Integrative and Comparative Biology__, 2025  
</details>

> 🎓 philentropy has been used in **dozens of peer-reviewed publications** to quantify distances, divergences, and similarities in complex biological and computational datasets.

---

## 🧠 Important Functions

### Distance Measures
- `distance()` – Implements 46 probability distance/similarity measures  
- `getDistMethods()` – Get all method names for `distance`  
- `dist.diversity()` – Computes distance diversity between PDFs  
- `estimate.probability()` – Estimate probability vectors from counts  

### Information Theory
- `H()` – Shannon’s Entropy `H(X)`  
- `JE()` – Joint Entropy `H(X,Y)`  
- `CE()` – Conditional Entropy `H(X|Y)`  
- `MI()` – Mutual Information `I(X,Y)`  
- `KL()` – Kullback–Leibler Divergence  
- `JSD()` – Jensen–Shannon Divergence  
- `gJSD()` – Generalized Jensen–Shannon Divergence  

---

## 🗞️ NEWS

Find the current status and version history in the  
👉 [NEWS section](https://drostlab.github.io/philentropy/news/index.html).

---

## 🧩 Appendix — full references

> - __A transcriptomic hourglass in brown algae__  
>   JS Lotharukpong, M Zheng, R Luthringer et al. – __Nature__, 2024
>
> - __Annelid functional genomics reveal the origins of bilaterian life cycles__  
>   FM Martín-Zamora, Y Liang, K Guynes et al. – __Nature__, 2023
>
> - __An atlas of gene regulatory elements in adult mouse cerebrum__  
>   YE Li, S Preissl, X Hou, Z Zhang, K Zhang et al. – __Nature__, 2021
>
> - __Convergent somatic mutations in metabolism genes in chronic liver disease__  
>   S Ng, F Rouhani, S Brunner, N Brzozowska et al. – __Nature__, 2021
>
> - __Antigen dominance hierarchies shape TCF1+ progenitor CD8 T cell phenotypes in tumors__  
>   ML Burger, AM Cruz, GE Crossland et al. – __Cell__, 2021
>
> - __High-content single-cell combinatorial indexing__  
>   R Mulqueen et al. – __Nature Biotechnology__, 2021
>
> - __A comparative atlas of single-cell chromatin accessibility in the human brain__  
>   YE Li, S Preissl, M Miller, ND Johnson, Z Wang et al. – __Science__, 2023
>
> - __Extinction at the end-Cretaceous and the origin of modern Neotropical rainforests__  
>   MR Carvalho, C Jaramillo et al. – __Science__, 2021
>
> - __sciCSR infers B cell state transition and predicts class-switch recombination dynamics using single-cell transcriptomic data__  
>   JCF Ng, G Montamat Garcia, AT Stewart et al. – __Nature Methods__, 2024
>
> - __HERMES: a molecular-formula-oriented method to target the metabolome__  
>   R Giné, J Capellades, JM Badia et al. – __Nature Methods__, 2021
>
> - __Epithelial zonation along the mouse and human small intestine defines five discrete metabolic domains__  
>   RK Zwick, P Kasparek, B Palikuqi et al. – __Nature Cell Biology__, 2024
>
> - __The genetic architecture of temperature adaptation is shaped by population ancestry and not by selection regime__  
>   KA Otte, V Nolte, F Mallard et al. – __Genome Biology__, 2021
>
> - __The Tug1 lncRNA locus is essential for male fertility__  
>   JP Lewandowski et al. – __Genome Biology__, 2020
>
> - __Decoding the gene regulatory network of endosperm differentiation in maize__  
>   Y Yuan, Q Huo, Z Zhang, Q Wang, J Wang et al. – __Nature Communications__, 2024
>
> - __A full-body transcription factor expression atlas with completely resolved cell identities in C. elegans__  
>   Y Li, S Chen, W Liu, D Zhao, Y Gao, S Hu, H Liu, Y Li et al. – __Nature Communications__, 2024
>
> - __Comprehensive mapping and modelling of the rice regulome landscape unveils the regulatory architecture underlying complex traits__  
>   T Zhu, C Xia, R Yu, X Zhou, X Xu, L Wang et al. – __Nature Communications__, 2024
>
> - __Transcriptional vulnerabilities of striatal neurons in human and rodent models of Huntington's disease__  
>   A Matsushima, SS Pineda, JR Crittenden et al. – __Nature Communications__, 2023
>
> - __Population structure in a fungal human pathogen is potentially linked to pathogenicity__  
>   EA Hatmaker, AE Barber, MT Drott et al. – __Nature Communications__, 2025
>
> - __Resolving the structure of phage–bacteria interactions in the context of natural diversity__  
>   KM Kauffman, WK Chang, JM Brown et al. – __Nature Communications__, 2022
>
> - __Gut microbiome-mediated metabolism effects on immunity in rural and urban African populations__  
>   M Stražar, GS Temba, H Vlamakis et al. – __Nature Communications__, 2021
>
> - __Aging, inflammation and DNA damage in the somatic testicular niche with idiopathic germ cell aplasia__  
>   M Alfano, AS Tascini, F Pederzoli et al. – __Nature Communications__, 2021
>
> - __Single cell census of human kidney organoids shows reproducibility and diminished off-target cells after transplantation__  
>   A Subramanian et al. – __Nature Communications__, 2019
>
> - __Pan-cancer human brain metastases atlas at single-cell resolution__  
>   X Xing, J Zhong, J Biermann, H Duan, X Zhang et al. – __Cancer Cell__, 2025
>
> - __The temporal progression of lung immune remodeling during breast cancer metastasis__  
>   CS McGinnis, Z Miao, D Superville, W Yao et al. – __Cancer Cell__, 2024
>
> - __Cross-tissue human fibroblast atlas reveals myofibroblast subtypes with distinct roles in immune modulation__  
>   Y Gao, J Li, W Cheng, T Diao, H Liu, Y Bo, C Liu et al. – __Cancer Cell__, 2024
>
> - __Gene module reconstruction identifies cellular differentiation processes and the regulatory logic of specialized secretion in zebrafish__  
>   Y Wang, J Liu, LY Du, JL Wyss, JA Farrell, AF Schier – __Developmental Cell__, 2025
>
> - __Large-scale chromatin reorganization reactivates placenta-specific genes that drive cellular aging__  
>   Z Liu, Q Ji, J Ren, P Yan, Z Wu, S Wang, L Sun, Z Wang et al. – __Developmental Cell__, 2022
>
> - __Integrated single-cell and spatial transcriptomic profiling reveals that CD177+ Tregs enhance immunosuppression through apoptosis and resistance to …__  
>   Y Liang, L Qiao, Q Qian, R Zhang, Y Li, X Xu et al. – __Oncogene__, 2025
>
> - __Conserved and unique features of terminal telomeric sequences in ALT-positive cancer cells__  
>   B Azeroglu, W Wu, R Pavani, RS Sandhu et al. – __eLife__, 2025
>
> - __Spotless, a reproducible pipeline for benchmarking cell type deconvolution in spatial transcriptomics__  
>   C Sang-Aram, R Browaeys, R Seurinck, Y Saeys – __eLife__, 2024
>
> - __Loss of adaptive capacity in asthmatic patients revealed by biomarker fluctuation dynamics after rhinovirus challenge__  
>   A Sinha et al. – __eLife__, 2019
>
> - __Staphylococci in high resolution: Capturing diversity within the human nasal microbiota__  
>   AC Ingham, DYK Ng, S Iversen, CM Liu et al. – __Cell Reports__, 2025
>
> - __Triple network dynamics and future alcohol consumption in adolescents__  
>   CC McIntyre, M Khodaei, RG Lyday et al. – __Alcohol: Clinical and Experimental Research__, 2025
>
> - __Benchmarking 13 tools for mutational signature attribution, including a new and improved algorithm__  
>   N Jiang, Y Wu, SG Rozen – __Briefings in Bioinformatics__, 2025
>
> - __Plant species as ecological engineers of microtopography in a temperate sedge-grass marsh__  
>   J Dušek, J Novotný, B Navrátilová et al. – __Scientific Reports__, 2025
>
> - __Resolution of MALDI-TOF compared to whole genome sequencing for identification of Bacillus species isolated from cleanrooms at NASA Johnson Space Center__  
>   F Mazhari, AB Regberg, CL Castro, MG LaMontagne – __Frontiers in Microbiology__, 2025
>
> - __Every Hue Has Its Fan Club: Diverse Patterns of Color-Dependent Flower Visitation across Lepidoptera__  
>   D Kutcherov, EL Westerman – __Integrative and Comparative Biology__, 2025
>
> - __An in vivo CRISPR screen in chick embryos reveals a role for MLLT3 in specification of neural cells from the caudal epiblast__  
>   ARG Libby, T Rito, A Radley, J Briscoe – __Development__, 2025
>
> - __Single-Cell Analyses Reveal a Functionally Heterogeneous Exhausted CD8+ T-cell Subpopulation That Is Correlated with Response to Checkpoint Therapy in …__  
>   KM Mahuron, O Shahid, P Sao, C Wu et al. – __Cancer Research__, 2025
>
> - __ETS1‐Driven Nucleolar Stress Orchestrates OLR1+ Macrophage Crosstalk to Sustain Immunosuppressive Microenvironment in Clear Cell Renal Cell Carcinoma__  
>   L Xiao, Z Zhang, T Li, Y Jiang, Y Liu, J Wang, W Tang – __Human Mutation__, 2025
>
> - __Association Between Ocular Microbiomes of Children and Their Siblings and Parents__  
>   X Ling, Y Zhang, CHT Bui, HN Chan, POS Tam et al. – __Investigative Ophthalmology & Visual Science__, 2025
>
> - __Benefits and challenges of host depletion methods in profiling the upper and lower respiratory microbiome__  
>   C Wang, L Zhang, C Kan, J He, W Liang, R Xia et al. – __Biofilms and Microbiomes__, 2025
>
> - __Unsettled Times: Music Discovery Reveals Divergent Cultural Responses to War__  
>   H Lee, M Anglada-Tort, O Sobchuk et al. – __PsyArXiv Preprints__, 2025
>
> - __q-Generalization of Nakagami distribution with applications__  
>   N Kumar, A Dixit, V Vijay – __Japanese Journal of Statistics and Data Science__, 2025
>
> - __The power of visualizing distributional differences: formal graphical n-sample tests__  
>   K Konstantinou, T Mrkvička, M Myllymäki – __Computational Statistics__, 2025
>
> - __Topic Modeling of Positive and Negative Reviews of Soulslike Video Games__  
>   T Guzsvinecz – __Computers__, 2025
>
> - __Basic Statistical Inference__  
>   M Nguyen – __Foundations of Data Analysis__, 2025