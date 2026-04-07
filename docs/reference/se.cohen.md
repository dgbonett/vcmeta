# Computes the standard error for Cohen's d

Computes the standard error of Cohen's d using only the two sample sizes
and an estimate of Cohen's d. Cohen's d and its standard error assume
equal variances. The estimate of Cohen's d, with the standard error
output from this function, can be used as input in the
[meta.ave.gen](https://dgbonett.github.io/vcmeta/reference/meta.ave.gen.md),
[meta.lc.gen](https://dgbonett.github.io/vcmeta/reference/meta.lc.gen.md),
and
[meta.lm.gen](https://dgbonett.github.io/vcmeta/reference/meta.lm.gen.md)
functions in applications where different types of compatible
standardized mean differences are used in the meta-analysis. If the
means, standard deviations, and sample sizes for the two groups are
available, use the
[se.stdmean2](https://dgbonett.github.io/vcmeta/reference/se.stdmean2.md)
function which does not assume equal variances. The standard error for
Cohen's d can be very inaccurate if the variances are unequal and the
sample sizes are unequal.

For more details, see Chapter 1 of Bonett (2021, Volume 5)

## Usage

``` r
se.cohen(d, n1, n2)
```

## Arguments

- d:

  estimated Cohen's d

- n1:

  sample size for group 1

- n2:

  sample size for group 2

## Value

Returns a one-row matrix:

- Estimate - Cohen's d (from input)

- SE - standard error

## References

Bonett DG (2021). *Statistical Methods for Psychologists, Vol 1-5,
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
se.cohen(.782, 35, 50)
#>             Estimate      SE
#> Cohen's d:     0.782 0.22887

# Should return: 
#            Estimate      SE
# Cohen's d:    0.782 0.22887

```
