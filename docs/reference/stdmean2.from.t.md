# Computes Cohen's d from pooled-variance t statistic

This function computes Cohen's d for a 2-group design (which is a
standardized mean difference with a weighted variance standardizer)
using a pooled-variance independent-samples t statistic and the two
sample sizes. This function also computes the standard error for Cohen's
d. The Cohen's d estimate and standard error assume equality of
population variances.

## Usage

``` r
stdmean2.from.t(t, n1, n2)
```

## Arguments

- t:

  pooled-variance t statistic

- n1:

  sample size for group 1

- n2:

  sample size for group 2

## Value

Returns Cohen's d and its equal-variance standard error

## Examples

``` r
stdmean2.from.t(3.27, 25, 25)
#>             Estimate     SE
#> Cohen's d:     0.925 0.2982

# Should return:
#           Estimate     SE
# Cohen's d    0.925 0.2982
```
