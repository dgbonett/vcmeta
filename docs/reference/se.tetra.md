# Computes the standard error for a tetrachoric correlation approximation

Computes an estimate of a tetrachoric correlation approximation and its
standard error using the frequency counts from a 2 x 2 contingency table
for two artifically dichotomous variables. A tetrachoric approximation
could be compatible with a Pearson correlation in a meta-analysis. The
tetrachoric approximation and the standard error from this function can
be used as input in the
[meta.ave.cor.gen](https://dgbonett.github.io/vcmeta/reference/meta.ave.cor.gen.md)
function in a meta-analysis where some studies have reported Pearson
correlations between quantitative variables x and y and other studies
have reported a 2 x 2 contingency table for dichotomous measurements of
variables x and y.

For more details, see Chapter 1 of Bonett (2021, Volume 5)

## Usage

``` r
se.tetra(f00, f01, f10, f11)
```

## Arguments

- f00:

  number of participants with y = 0 and x = 0

- f01:

  number of participants with y = 0 and x = 1

- f10:

  number of participants with y = 1 and x = 0

- f11:

  number of participants with y = 1 and x = 1

## Value

Returns a 1-row matrix. The columns are:

- Estimate - estimated tetrachoric approximation

- SE - standard error

## References

Bonett DG, Price RM (2005). “Inferential methods for the tetrachoric
correlation coefficient.” *Journal of Educational and Behavioral
Statistics*, **30**(2), 213–225. ISSN 1076-9986,
[doi:10.3102/10769986030002213](https://doi.org/10.3102/10769986030002213)
.

Bonett DG (2021). *Statistical Methods for Psychologists, Vol 1-5,
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
se.tetra(46, 15, 54, 85)
#>               Estimate     SE
#> Tetrachoric:     0.514 0.0936

# Should return:
#               Estimate     SE 
# Tetrachoric:     0.514 0.0936

```
