# Computes the standard error for a Spearman correlation

Computes the Bonett-Wright standard error of a Spearman correlation
using the estimated Spearman correlation and sample size. The standard
error from this function can be used as input in the
[meta.ave.cor.gen](https://dgbonett.github.io/vcmeta/reference/meta.ave.cor.gen.md)
function in applications where a combination of different types of
compatible correlations are used in the meta-analysis.

For more details, see Chapter 1 of Bonett (2021, Volume 5)

## Usage

``` r
se.spear(cor, n)
```

## Arguments

- cor:

  estimated Spearman correlation

- n:

  sample size

## Value

Returns a one-row matrix:

- Estimate - Spearman correlation (from input)

- SE - standard error

## References

Bonett DG, Wright TA (2000). “Sample size requirements for estimating
Pearson, Kendall and Spearman correlations.” *Psychometrika*, **65**(1),
23–28. ISSN 0033-3123,
[doi:10.1007/BF02294183](https://doi.org/10.1007/BF02294183) .

Bonett DG (2021). *Statistical Methods for Psychologists, Vol 1-5,
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
se.spear(.427, 55)
#>                        Estimate      SE
#> Spearman correlation:     0.427 0.11845

# Should return: 
#                       Estimate      SE
# Spearman correlation:    0.427 0.11845

```
