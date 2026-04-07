# Confidence interval for a linear contrast of means

Computes the estimate, standard error, and confidence interval for a
linear contrast of means from two or more studies. This function will
use either an unequal variance (recommended) or an equal variance
method. A Satterthwaite adjustment to the degrees of freedom is used
with the unequal variance method.

For more details, see Setion 3.2 of Bonett (2021, Volume 5).

## Usage

``` r
meta.lc.mean1(alpha, m, sd, n, v, eqvar = FALSE)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- m:

  vector of estimated means

- sd:

  vector of estimated standard deviations

- n:

  vector of sample sizes

- v:

  vector of contrast coefficients

- eqvar:

  - FALSE for unequal variance method

  - TRUE for equal variance method

## Value

Returns 1-row matrix with the following columns:

- Estimate - estimated linear contrast

- SE - standard error

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

- df - degrees of freedom

## References

Bonett DG (2021). *Statistical Methods for Psychologists, Vol 1-5,
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
m <- c(33.5, 37.9, 38.0, 44.1)
sd <- c(3.84, 3.84, 3.65, 4.98)
n <- c(10, 10, 10, 10)
v <- c(.5, .5, -.5, -.5)
meta.lc.mean1(.05, m, sd, n, v, eqvar = FALSE)
#>          Estimate       SE        LL        UL    df
#> Contrast    -5.35 1.300136 -7.993588 -2.706412 33.52

# Should return:
#          Estimate       SE        LL        UL    df
# Contrast    -5.35 1.300136 -7.993583 -2.706417 33.52

```
