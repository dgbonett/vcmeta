# Computes the standard error for the average of two Pearson correlations with no variables in common that have been estimated from the same sample

In a study that reports the sample size and six correlations (cor12,
cor34, cor13, cor14, cor23, and cor24) where variables 1 and 3 are
different measurements of one attribute and variables 2 and 4 are
different measurements of a second attribute, this function can be used
to compute the average of cor12 and cor34 and its standard error. Note
that cor12 and cor34 have no variable in common (i.e., no "overlapping"
variable). The average correlation and the standard error from this
function can be used as input in the
[meta.ave.cor.gen](https://dgbonett.github.io/vcmeta/reference/meta.ave.cor.gen.md)
function in a meta-analysis where some studies have reported cor12 and
other studies have reported cor34.

For more details, see Chapter 1 of Bonett (2021, Volume 5)

## Usage

``` r
se.ave.cor.nonover(cor12, cor34, cor13, cor14, cor23, cor24, n)
```

## Arguments

- cor12:

  estimated correlation between variables 1 and 2

- cor34:

  estimated correlation between variables 3 and 4

- cor13:

  estimated correlation between variables 1 and 3

- cor14:

  estimated correlation between variables 1 and 4

- cor23:

  estimated correlation between variables 2 and 3

- cor24:

  estimated correlation between variables 2 and 4

- n:

  sample size

## Value

Returns a two-row matrix. The first row gives results for the average
correlation and the second row gives the results with a Fisher
transformation. The columns are:

- Estimate - estimated average of cor12 and cor34

- SE - standard error

- VAR(cor12) - variance of cor12

- VAR(cor34) - variance of cor34

- COV(cor12,cor34) - covariance of cor12 and cor34

## References

Bonett DG (2021). *Statistical Methods for Psychologists, Vol 1-5,
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
se.ave.cor.nonover(.357, .398, .755, .331, .347, .821, 100)
#>               Estimate      SE VAR(cor12)  VAR(cor34) COV(cor12,cor34)
#> Correlation:    0.3775 0.07769 0.00784892 0.007301895      0.004495714
#> Fisher:         0.3971 0.09060 0.01030928 0.010309278      0.006122153

# Should return:
#              Estimate      SE VAR(cor12)  VAR(cor34) COV(cor12,cor34)
# Correlation:   0.3775 0.07769 0.00784892 0.007301895      0.004495714
# Fisher:        0.3971 0.09060 0.01030928 0.010309278      0.006122153

```
