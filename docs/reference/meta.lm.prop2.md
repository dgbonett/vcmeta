# Meta-regression analysis for 2-group proportion differences

This function estimates the intercept and slope coefficients in a
meta-regression model where the dependent variable is a 2-group
proportion difference. The estimates are OLS estimates with robust
standard errors that accommodate residual heteroscedasticity.

For more details, see Section 3.4 of Bonett (2021, Volume 5).

## Usage

``` r
meta.lm.prop2(alpha, f1, f2, n1, n2, X)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- f1:

  vector of group 1 frequency counts

- f2:

  vector of group 2 frequency counts

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

- z - z-value

- p - p-value

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## References

Bonett DG, Price RM (2014). “Meta-analysis methods for risk
differences.” *British Journal of Mathematical and Statistical
Psychology*, **67**(3), 371–387. ISSN 00071102,
[doi:10.1111/bmsp.12024](https://doi.org/10.1111/bmsp.12024) .

Bonett DG (2021). *Statistical Methods for Psychologists, Vol 1-5,
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
f1 <- c(24, 40, 93, 14, 5)
f2 <- c(12, 9, 28, 3, 1)
n1 <- c(204, 201, 932, 130, 77)
n2 <- c(106, 103, 415, 132, 83)
x1 <- c(4, 4, 5, 3, 26)
x2 <- c(1, 1, 1, 0, 0)
X <- matrix(cbind(x1, x2), 5, 2)
meta.lm.prop2(.05, f1, f2, n1, n2, X)
#>        Estimate          SE      z     p          LL          UL
#> b0  0.089756283 0.034538077  2.599 0.009  0.02206290 0.157449671
#> b1 -0.001447968 0.001893097 -0.765 0.444 -0.00515837 0.002262434
#> b2 -0.034670988 0.034125708 -1.016 0.310 -0.10155615 0.032214170

# Should return:
#        Estimate          SE      z     p          LL          UL
# b0  0.089756283 0.034538077  2.599 0.009  0.02206290 0.157449671
# b1 -0.001447968 0.001893097 -0.765 0.444 -0.00515837 0.002262434
# b2 -0.034670988 0.034125708 -1.016 0.310 -0.10155615 0.032214170

```
