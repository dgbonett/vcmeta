# Computes the estimate and standard error for a paired-samples proportion difference

This function computes a standard error of a paired-samples proportion
difference using the frequency counts from a 2 x 2 contingency table and
the number of studies in the planned meta-analysis. The effect size
estimate and standard error output from this function can be used as
input in the
[meta.ave.gen](https://dgbonett.github.io/vcmeta/reference/meta.ave.gen.md),
[meta.lc.gen](https://dgbonett.github.io/vcmeta/reference/meta.lc.gen.md),
and
[meta.lm.gen](https://dgbonett.github.io/vcmeta/reference/meta.lm.gen.md)
functions in applications where compatible proportion differences from a
combination of 2-group and paired-samples studies are used in the
meta-analysis.

For more details, see Chapter 1 of Bonett (2021, Volume 5)

## Usage

``` r
se.prop.ps(f00, f01, f10, f11, m)
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

- m:

  number of studies in planned meta-analysis

## Value

Returns a one-row matrix:

- Estimate - adjusted estimate of proportion difference

- SE - standard error

## References

Bonett DG, Price RM (2012). “Adjusted wald confidence interval for a
difference of binomial proportions based on paired data.” *Journal of
Educational and Behavioral Statistics*, **37**(4), 479–488. ISSN
1076-9986,
[doi:10.3102/1076998611411915](https://doi.org/10.3102/1076998611411915)
.

Bonett DG (2021). *Statistical Methods for Psychologists, Vol 1-5,
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
se.prop.ps(16, 64, 5, 15, 4)
#>                          Estimate        SE
#> Proportion difference:  0.5870647 0.0587513

# Should return:
#                          Estimate        SE
# Proportion difference:  0.5870647 0.0587513

```
