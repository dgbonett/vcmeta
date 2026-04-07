# Meta-regression analysis for paired-samples mean differences

This function estimates the intercept and slope coefficients in a
meta-regression model where the dependent variable is a paired-samples
mean difference. The estimates are OLS estimates with robust standard
errors that accommodate residual heteroscedasticity.

For more details, see Section 3.4 of Bonett (2021, Volume 5).

## Usage

``` r
meta.lm.mean.ps(alpha, m1, m2, sd1, sd2, cor, n, X)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- m1:

  vector of estimated means for measurement 1

- m2:

  vector of estimated means for measurement 2

- sd1:

  vector of estimated SDs for measurement 1

- sd2:

  vector of estimated SDs for measurement 2

- cor:

  vector of estimated correlations for paired measurements

- n:

  vector of sample sizes

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
n <- c(65, 30, 29, 45, 50)
cor <- c(.87, .92, .85, .90, .88)
m1 <- c(20.1, 20.5, 19.3, 21.5, 19.4)
m2 <- c(10.4, 10.2, 8.5, 10.3, 7.8)
sd1 <- c(9.3, 9.9, 10.1, 10.5, 9.8)
sd2 <- c(7.8, 8.0, 8.4, 8.1, 8.7)
x1 <- c(2, 3, 3, 4, 4)
X <- matrix(x1, 5, 1)
meta.lm.mean.ps(.05, m1, m2, sd1, sd2, cor, n, X)
#>    Estimate        SE     t     p        LL        UL  df
#> b0     8.00 1.2491990 6.404 0.000 5.5378833 10.462117 217
#> b1     0.85 0.3796019 2.239 0.026 0.1018213  1.598179 217

# Should return:
#    Estimate        SE     t     p        LL        UL  df
# b0     8.00 1.2491990 6.404 0.000 5.5378833 10.462117 217
# b1     0.85 0.3796019 2.239 0.026 0.1018213  1.598179 217

```
