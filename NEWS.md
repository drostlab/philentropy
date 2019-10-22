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
