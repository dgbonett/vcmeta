# Computes the standard error for the average of two Pearson correlations with one variable in common that have been estimated from the same sample

In a study that reports the sample size and three correlations (cor12,
cor13, and cor23 where variable 1 is called the "overlapping" variable),
and variables 2 and 3 are different measurements of the same attribute,
this function can be used to compute the average of cor12 and cor13 and
its standard error. The average correlation and the standard error from
this function can be used as input in the
[meta.ave.cor.gen](https://dgbonett.github.io/vcmeta/reference/meta.ave.cor.gen.md)
function in a meta-analysis where some studies have reported cor12 and
other studies have reported cor13.

For more details, see Chapter 1 of Bonett (2021, Volume 5)

## Usage

``` r
se.ave.cor.over(cor12, cor13, cor23, n)
```

## Arguments

- cor12:

  estimated correlation between variables 1 and 2

- cor13:

  estimated correlation between variables 1 and 3

- cor23:

  estimated correlation between variables 2 and 3

- n:

  sample size

## Value

Returns a two-row matrix. The first row gives results for the average
correlation and the second row gives the results with a Fisher
transformation. The columns are:

- Estimate - estimated average of cor12 and cor13

- SE - standard error

- VAR(cor12) - variance of cor12

- VAR(cor13) - variance of cor13

- COV(cor12,cor13) - covariance of cor12 and cor13

## References

Bonett DG (2021). *Statistical Methods for Psychologists, Vol 1-5,
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
se.ave.cor.over(.462, .518, .755, 100)
#>               Estimate      SE  VAR(cor12) VAR(cor13) COV(cor12,cor13)
#> Correlation:    0.4900 0.07087 0.006378045 0.00551907      0.004097553
#> Fisher:         0.5361 0.09327 0.010309278 0.01030928      0.007119936

# Should return:
#              Estimate      SE  VAR(cor12) VAR(cor13) COV(cor12,cor13)
# Correlation:   0.4900 0.07087 0.006378045 0.00551907      0.004097553
# Fisher:        0.5361 0.09327 0.010309278 0.01030928      0.007119936

```
