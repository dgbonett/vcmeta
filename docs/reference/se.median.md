# Computes the standard error for a median

Computes the standard error of a median using a confidence interval for
the distribution-free interval (see Snedecor & Cochran, 1980), this
function computes a Price-Bonett standard error which is known to be
very accurate. If any other type of confidence interval is used, this
function computes an approximate standard error.

In a 2-group design, this function can be used to compute the standard
error of the median in each group. Then the difference in estimated
medians in the two groups can be used as the effect size and the
standard error of the difference can be set to the square root of the
sum of the squared standard errors from each group. Single-group medians
or median differences and their standard errors can be used as input in
the
[meta.ave.gen](https://dgbonett.github.io/vcmeta/reference/meta.ave.gen.md),
[meta.lc.gen](https://dgbonett.github.io/vcmeta/reference/meta.lc.gen.md),
and and
[meta.lm.gen](https://dgbonett.github.io/vcmeta/reference/meta.lm.gen.md)
functions.

## Usage

``` r
se.median(alpha, LL, UL, n, type)
```

## Arguments

- alpha:

  alpha value for 1-alpha confidence interval

- LL:

  lower limit of confidence interval

- UL:

  upper limit of confidence interval

- n:

  sample size

- type:

  - set to 1 for classical confidence interval

  - set to 2 for any other type of confidence interval

## Value

Returns the estimated standard error of the sample median

## References

Price RM, Bonett DG (2001). “Estimating the variance of the sample
median.” *Journal of Computation and Simulation*, **68**(3), 295–305.
[doi:10.1080/00949650108812071](https://doi.org/10.1080/00949650108812071)
.

Snedecor GW, Cochran WG (1980). *Statistical methods*, 7th edition. ISU
University Pres, Ames, Iowa.

## Examples

``` r
se.median(.05, 47.21, 68.68, 35, 1)
#>        SE
#>  5.252114

# Should return:
#        SE
#  5.252114

```
