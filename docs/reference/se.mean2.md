# Computes the standard error for a 2-group mean difference

Computes the standard error of a 2-group mean difference using the
estimated means, estimated standard deviations, and sample sizes. The
effect size estimate and standard error output from this function can be
used as input in the
[meta.ave.gen](https://dgbonett.github.io/vcmeta/reference/meta.ave.gen.md),
[meta.lc.gen](https://dgbonett.github.io/vcmeta/reference/meta.lc.gen.md),
and
[meta.lm.gen](https://dgbonett.github.io/vcmeta/reference/meta.lm.gen.md)
functions in applications where compatible mean differences from a
combination of 2-group and paired-samples experiments are used in the
meta-analysis. Equality of variances is not asumed.

For more details, see Chapter 1 of Bonett (2021, Volume 5)

## Usage

``` r
se.mean2(m1, m2, sd1, sd2, n1, n2)
```

## Arguments

- m1:

  estimated mean for group 1

- m2:

  estimated mean for group 2

- sd1:

  estimated standard deviation for group 1

- sd2:

  estimated standard deviation for group 2

- n1:

  sample size for group 1

- n2:

  sample size for group 2

## Value

Returns a one-row matrix:

- Estimate - estimated mean difference

- SE - standard error

## References

Bonett DG (2021). *Statistical Methods for Psychologists, Vol 1-5,
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
se.mean2(21.93, 16.11, 3.82, 3.21, 40, 40)
#>                   Estimate        SE
#> Mean difference:      5.82 0.7889312

# Should return:
#                   Estimate        SE
# Mean difference:      5.82 0.7889312

```
