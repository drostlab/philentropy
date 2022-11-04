## Version 0.7.0

### New Features

### Updates

- the [Distances](https://drostlab.github.io/philentropy/articles/Distances.html) vignette now has a fixed documentation for the benchmarking of low-level distance functions. Many thanks to (@Nowosad) #30 
- in `../src/correlation.h` adjustment of use of logical operators rather than Wbitwise (`| -> or`) which otherwises raises warnings in `clang14`
- vector element limit is now extended to long vectors for all distance measures by declaring `R_xlen_t` instead of `int` during indexing.

## Version 0.6.0

### New Features

- `distance()` and all other individual information theory functions
receive a new argument `epsilon` with default value `epsilon = 0.00001` to treat cases where in individual distance or similarity computations 
yield `x / 0` or `0 / 0`. Instead of a hard coded epsilon, users can now set `epsilon` according to their input vectors. (Many thanks to Joshua McNeill #26 for this great question). 
- three new functions `dist_one_one()`, `dist_one_many()`, `dist_many_many()` are added. They are fairly flexible intermediaries between `distance()` and single distance functions. `dist_one_one()` expects two vectors (probability density functions) and returns a single value. `dist_one_many()` expects one vector (a probability density function) and one matrix (a set of probability density functions), and returns a vector of values. `dist_many_many()` expects two matrices (two sets of probability density functions), and returns a matrix of values. (Many thanks to 
Jakub Nowosad, see #27, #28, and [New Vignette Many_Distance](https://drostlab.github.io/philentropy/articles/Many_Distances.html))

### Updates

- a new Vignette [Comparing many probability density functions](https://drostlab.github.io/philentropy/articles/Many_Distances.html) (Many thanks to 
Jakub Nowosad)
- `dplyr` package dependency was removed and replaced by the `poorman`
due to the heavy dependency burden of `dplyr`, since `philentropy`
only used `dplyr::between()` which is now `poorman::between()` (Many thanks to Patrice Kiener for this suggestion)
- `distance(..., as.dist.obj = TRUE)` now returns the same values as `stats::dist()` when working with 2 dimensional input matrices (2 vector inputs) (see #29) (Many thanks to 
Jakub Nowosad (@Nowosad))
Example:

```r
library(philentropy)

m1 = matrix(c(1, 2), ncol = 1)

dist(m1)
#> 1
#> 2 1
distance(m1, as.dist.obj = TRUE)
#> Metric: 'euclidean'; comparing: 2 vectors.
#> 1
#> 2 1
```


## Version 0.5.0

### New Features

- the `distance()` function receives a new argument `mute.message` allowing users to mute  message printing when running large-scale distance computations.
Example:

```r
distance(rbind(1:10/sum(1:10), 20:29/sum(20:29)), 
         method = "euclidean", 
         mute.message = TRUE)
```

- adding `markdown` dependency to `DESCRIPTION` ([find details here](https://github.com/yihui/knitr/issues/1864))

## Version 0.4.0

### New Features

- the `distance()` function receives a new argument `use.row.names` to enable passing the row names from the input probability or count matrix to the output distance matrix

- the `distance()` function can now handle `data.table` and `tibble` input #16

- adding new functionality and arguments `as.dist.obj`, `diag`, and `upper` to `philentropy::distance()` to allow users to retrieve a `stats::dist()` object when working with `philentropy::distance()` (Many thanks to Hugo Tavares #18 - see also #13)
When using `philentropy::distance(..., as.dist.obj = TRUE)` users can now directly pass the `distance()` output into `hclust`:

Before:
```r
ProbMatrix <- rbind(1:10/sum(1:10), 20:29/sum(20:29),30:39/sum(30:39))
dist.mat <- distance(ProbMatrix, method = "jaccard")
true.dist.mat <- as.dist(dist.mat)
clust.res <- hclust(true.dist.mat, method = "complete")
clust.res
```

```
Call:
hclust(d = true.dist.mat, method = "complete")

Cluster method   : complete 
Number of objects: 3 
```

Now:

```r
ProbMatrix <- rbind(1:10/sum(1:10), 20:29/sum(20:29),30:39/sum(30:39))
dist.mat <- distance(ProbMatrix, method = "jaccard", as.dist.obj = TRUE)
clust.res <- hclust(true.dist.mat, method = "complete")
clust.res
```

```
Call:
hclust(d = true.dist.mat, method = "complete")

Cluster method   : complete 
Number of objects: 3 
```

### Bug fixes

- fixing a bug in `gJSD()` which tested transposed matrix rows rather than transposed matrix columns for sum > 1 (see issue #17 ; many thanks to @wkc1986)

## Version 0.3.0

### New functionality
- exporting all Rcpp distance measure functions individually (see issue #9), this
enables access to much faster computations (see micro benchmarks at https://hajkd.github.io/philentropy/articles/Distances.html)

### Bug fixes

- fixing bug which caused that KL distance returns NaN when P == 0 (see issue #10; Many thanks to @KaiserDominici)

- fixing bug which caused stack overflow when computing distance matrices with many rows (see issue #7; Many thanks to @wkc1986 and @elbamos)

- fixing bug in `gJSD()` where an `rbind()` input matrix is not properly transposed (Many thanks to @vrodriguezf; see issue #14) 


### New Features

- `gJSD()` receives new argument `est.prob` to enable empirical estimation of probability vectors from input count vectors (non-probabilistic vectors) 

- Jaccard and Tanimoto similarity measures now return `0` instead of `NAN` when probability vectors contain zeros (Many thanks to @JonasMandel; see issue #15)


## Version 0.2.0

### Bug fixes
- Fixing bug that caused `jensen-shannon` computations to compute wrong values when `0 values` were present in the input vectors (see issue #4 ; Many thanks to @wkc1986)
- Fixing bug that caused `jensen-difference` computations to compute wrong values when `0 values` were present in the input vectors
- Fixing bugs in all distance metrics when handing 0/0, 0/x or x/0 cases

## Version 0.1.0

### New Features

- new message system
- extending documentation

### Bug fixes

- Fixing bug that caused that `JSD()` gives NaN when any probability is 0 - see https://github.com/HajkD/philentropy/issues/1 (Thanks to William Kurtis Chang)

## Version 0.0.2

### Bug fixes

- Fixing C++ memory leaks in `dist.diversity()` and `distance()` when check for `colSums(x) > 1.001` was peformed (leak was found with `rhub::check_with_valgrind()`)

## Version 0.0.1

Initial submission version.
