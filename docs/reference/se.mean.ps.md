# Computes the standard error for a paired-samples mean difference

Computes the standard error of a paired-samples mean difference using
the estimated means, estimated standard deviations, estimated Pearson
correlation, and sample size. The effect size estimate and standard
error output from this function can be used as input in the
[meta.ave.gen](https://dgbonett.github.io/vcmeta/reference/meta.ave.gen.md),
[meta.lc.gen](https://dgbonett.github.io/vcmeta/reference/meta.lc.gen.md),
and
[meta.lm.gen](https://dgbonett.github.io/vcmeta/reference/meta.lm.gen.md)
functions in applications where compatible mean differences from a
combination of 2-group and paired-samples experiments are used in the
meta-analysis. Equality of variances is not assumed.

For more details, see Chapter 1 of Bonett (2021, Volume 5)

## Usage

``` r
se.mean.ps(m1, m2, sd1, sd2, cor, n)
```

## Arguments

- m1:

  estimated mean for measurement 1

- m2:

  estimated mean for measurement 2

- sd1:

  estimated standard deviation for measurement 1

- sd2:

  estimated standard deviation for measurement 2

- cor:

  estimated correlation for measurements 1 and 2

- n:

  sample size

## Value

Returns a one-row matrix:

- Estimate - estimated mean difference

- SE - standard error

## References

Bonett DG (2021). *Statistical Methods for Psychologists, Vol 1-5,
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
se.mean.ps(23.94, 25.12, 1.76, 2.01, .78, 25)
#>                   Estimate        SE
#> Mean difference:     -1.18 0.2544833

# Should return:
#                   Estimate        SE
# Mean difference:     -1.18 0.2544833

```
