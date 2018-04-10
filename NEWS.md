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
