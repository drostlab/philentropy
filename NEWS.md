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
