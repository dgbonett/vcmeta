# Meta-regression analysis for 2-group mean differences

This function estimates the intercept and slope coefficients in a
meta-regression model where the dependent variable is a 2-group mean
difference. The estimates are OLS estimates with robust standard errors
that accommodate residual heteroscedasticity.

For more details, see Section 3.4 of Bonett (2021, Volume 5).

## Usage

``` r
meta.lm.mean2(alpha, m1, m2, sd1, sd2, n1, n2, X)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- m1:

  vector of estimated means for group 1

- m2:

  vector of estimated means for group 2

- sd1:

  vector of estimated SDs for group 1

- sd2:

  vector of estimated SDs for group 2

- n1:

  vector of group 1 sample sizes

- n2:

  vector of group 2 sample sizes

- X:

  matrix of predictor values

## Value

Returns a matrix. The first row is for the intercept with one additional
row per predictor. The matrix has the following columns:

- Estimate - OLS estimate

- SE - standard error

- t - t-value

- p - p-value

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

- df - degrees of freedom

## References

Bonett DG (2009). “Meta-analytic interval estimation for standardized
and unstandardized mean differences.” *Psychological Methods*,
**14**(3), 225–238. ISSN 1939-1463,
[doi:10.1037/a0016619](https://doi.org/10.1037/a0016619) .

Bonett DG (2021). *Statistical Methods for Psychologists, Vol 1-5,
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
n1 <- c(65, 30, 29, 45, 50)
n2 <- c(67, 32, 31, 20, 52)
m1 <- c(31.1, 32.3, 31.9, 29.7, 33.0)
m2 <- c(34.1, 33.2, 30.6, 28.7, 26.5)
sd1 <- c(7.1, 8.1, 7.8, 6.8, 7.6)
sd2 <- c(7.8, 7.3, 7.5, 7.2, 6.8)
x1 <- c(4, 6, 7, 7, 8)
x2 <- c(1, 0, 0, 0, 1)
X <- matrix(cbind(x1, x2), 5, 2)
meta.lm.mean2(.05, m1, m2, sd1, sd2, n1, n2, X)
#>    Estimate        SE      t     p         LL        UL  df
#> b0   -15.20 3.4097610 -4.458 0.000 -21.902415 -8.497585 418
#> b1     2.35 0.4821523  4.874 0.000   1.402255  3.297745 418
#> b2     2.85 1.5358109  1.856 0.064  -0.168875  5.868875 418

# Should return:
#    Estimate        SE      t     p         LL        UL  df
# b0   -15.20 3.4097610 -4.458 0.000 -21.902415 -8.497585 418
# b1     2.35 0.4821523  4.874 0.000   1.402255  3.297745 418
# b2     2.85 1.5358109  1.856 0.064  -0.168875  5.868875 418

```
