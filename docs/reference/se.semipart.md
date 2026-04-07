# Computes the standard error for a semipartial correlation

Computes the standard error of a semipartial correlation using the
estimated semipartial correlation, sample size, and squared multiple
correlation for the full model. The full model includes the independent
variable of interest and all control variables. The effect size estimate
and standard error output from this function can be used as input in the
[meta.ave.cor.gen](https://dgbonett.github.io/vcmeta/reference/meta.ave.cor.gen.md)
function in applications where a combination of different types of
compatible correlations are used in the meta-analysis.

For more details, see Chapter 1 of Bonett (2021, Volume 5)

## Usage

``` r
se.semipart(cor, r2, n)
```

## Arguments

- cor:

  estimated semipartial correlation

- r2:

  estimated squared multiple correlation for a model that includes the
  IV and all control variables

- n:

  sample size

## Value

Returns a one-row matrix:

- Estimate - semipartial correlation (from input)

- SE - standard error

## References

Bonett DG (2021). *Statistical Methods for Psychologists, Vol 1-5,
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
se.semipart(.454, .25, 60)
#>                           Estimate      SE
#> Semipartial correlation:     0.454 0.10298

# Should return: 
#                           Estimate      SE
# Semipartial correlation:     0.454 0.10298

```
