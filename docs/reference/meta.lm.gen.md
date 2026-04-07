# Meta-regression analysis for any type of effect size

This function estimates the intercept and slope coefficients in a
meta-regression model where the dependent variable is any type of effect
size. The estimates are OLS estimates with robust standard errors that
accomodate residual heteroscedasticity.

For more details, see Section 3.4 of Bonett (2021, Volume 5).

## Usage

``` r
meta.lm.gen(alpha, est, se, X)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- est:

  vector of parameter estimates

- se:

  vector of standard errors

- X:

  matrix of predictor values

## Value

Returns a matrix. The first row is for the intercept with one additional
row per predictor. The matrix has the following columns:

- Estimate - OLS estimate

- SE - standard error

- z - z-value

- p - p-value

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## References

Bonett DG (2021). *Statistical Methods for Psychologists, Vol 1-5,
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
est <- c(4.1, 4.7, 4.9, 5.7, 6.6, 7.3)
se <- c(1.2, 1.5, 1.3, 1.8, 2.0, 2.6)
x1 <- c(10, 20, 30, 40, 50, 60)
x2 <- c(1, 1, 1, 0, 0, 0)
X <- matrix(cbind(x1, x2), 6, 2)
meta.lm.gen(.05, est, se, X)
#>      Estimate         SE      z     p         LL         UL
#> b0  3.5333333 4.37468253  0.808 0.419 -5.0408869 12.1075535
#> b1  0.0600000 0.09058835  0.662 0.508 -0.1175499  0.2375499
#> b2 -0.1666667 2.81139793 -0.059 0.953 -5.6769054  5.3435720

# Should return:
#      Estimate         SE      z     p         LL         UL
# b0  3.5333333 4.37468253  0.808 0.419 -5.0408869 12.1075535
# b1  0.0600000 0.09058835  0.662 0.508 -0.1175499  0.2375499
# b2 -0.1666667 2.81139793 -0.059 0.953 -5.6769054  5.3435720

```
