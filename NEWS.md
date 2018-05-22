## Version 0.2.0

### Bug fixes
- Fixing bug that caused `jensen-shannon` computations to compute wrong values when `0 values` were present in the input vectors (see issue #4 ; Many thanks to @wkc1986)
- Fixing bug that caused `jensen-difference` computations to compute wrong values when `0 values` were present in the input vectors
- Fixing bugs in all distance metrics when handing 0/0, 0/x or x/0 cases

### New Features


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
