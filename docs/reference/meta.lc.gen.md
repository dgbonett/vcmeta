# Confidence interval for a linear contrast of effect sizes

Computes the estimate, standard error, and confidence interval for a
linear contrast of any type of effect size from two or more studies.

For more details, see Setion 3.2 of Bonett (2021, Volume 5).

## Usage

``` r
meta.lc.gen(alpha, est, se, v)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- est:

  vector of parameter estimates

- se:

  vector of standard errors

- v:

  vector of contrast coefficients

## Value

Returns 1-row matrix with the following columns:

- Estimate - estimated linear contrast

- SE - standard error

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## References

Bonett DG (2021). *Statistical Methods for Psychologists, Vol 1-5,
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
est <- c(.55, .59, .44, .48, .26, .19)
se <- c(.054, .098, .029, .084, .104, .065)
v <- c(.5, .5, -.25, -.25, -.25, -.25)
meta.lc.gen(.05, est, se, v)
#>          Estimate         SE        LL        UL
#> Contrast   0.2275 0.06755461 0.0950954 0.3599046

# Should return: 
#          Estimate         SE        LL        UL
# Contrast   0.2275 0.06755461 0.0950954 0.3599046

```
