# Computes the standard error for a Pearson or partial correlation

Computes the standard error of a Pearson or partial correlation using
the estimated correlation, sample size, and number of control variables.
The correlation, along with the standard error output from this
function, can be used as input in the
[meta.ave.cor.gen](https://dgbonett.github.io/vcmeta/reference/meta.ave.cor.gen.md)
function in applications where a combination of different types of
compatible correlations are used in the meta-analysis.

For more details, see Chapter 1 of Bonett (2021, Volume 5)

## Usage

``` r
se.cor(cor, s, n)
```

## Arguments

- cor:

  estimated Pearson or partial correlation

- s:

  number of control variables (set to 0 for Pearson)

- n:

  sample size

## Value

Returns a one-row matrix:

- Estimate - Pearson or partial correlation (from input)

- SE - standard error

## References

Bonett DG (2008). “Meta-analytic interval estimation for bivariate
correlations.” *Psychological Methods*, **13**(3), 173–181. ISSN
1939-1463, [doi:10.1037/a0012868](https://doi.org/10.1037/a0012868) .

Bonett DG (2021). *Statistical Methods for Psychologists, Vol 1-5,
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
se.cor(.427, 0, 55)
#>               Estimate      SE
#> Correlation:     0.427 0.11339

# Should return: 
#               Estimate      SE
# Correlation:     0.427 0.11339

se.cor(.283, 4, 80)
#>               Estimate      SE
#> Correlation:     0.283 0.10767
# Should return: 
#               Estimate      SE
# Correlation:     0.283 0.10767

```
