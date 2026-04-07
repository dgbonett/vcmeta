# Meta-regression analysis for 2-group proportion ratios

This function estimates the intercept and slope coefficients in a
meta-regression model where the dependent variable is a 2-group log
proportion ratio. The estimates are OLS estimates with robust standard
errors that accommodate residual heteroscedasticity. The exponentiated
slope estimate for a predictor variable describes a multiplicative
change in the proportion ratio associated with a 1-unit increase in that
predictor variable, controlling for all other predictor variables in the
model.

For more details, see Section 3.4 of Bonett (2021, Volume 5).

## Usage

``` r
meta.lm.propratio2(alpha, f1, f2, n1, n2, X)
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

- exp(Estimate) - the exponentiated estimate

- exp(LL) - lower limit of the exponentiated confidence interval

- exp(UL) - upper limit of the exponentiated confidence interval

## References

Price RM, Bonett DG (2008). “Confidence intervals for a ratio of two
independent binomial proportions.” *Statistics in Medicine*, **27**(26),
5497–5508. ISSN 02776715,
[doi:10.1002/sim.3376](https://doi.org/10.1002/sim.3376) .

Bonett DG (2021). *Statistical Methods for Psychologists, Vol 1-5,
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
n1 <- c(204, 201, 932, 130, 77)
n2 <- c(106, 103, 415, 132, 83)
f1 <- c(24, 40, 93, 14, 5)
f2 <- c(12, 9, 28, 3, 1)
x1 <- c(4, 4, 5, 3, 26)
x2 <- c(1, 1, 1, 0, 0)
X <- matrix(cbind(x1, x2), 5, 2)
meta.lm.propratio2(.05, f1, f2, n1, n2, X)
#>         Estimate         SE      z     p          LL         UL exp(Estimate)
#> b0  1.4924887636 0.69172794  2.158 0.031  0.13672691 2.84825062     4.4481522
#> b1  0.0005759509 0.04999884  0.012 0.990 -0.09741998 0.09857188     1.0005761
#> b2 -1.0837844594 0.59448206 -1.823 0.068 -2.24894789 0.08137897     0.3383128
#>      exp(LL)   exp(UL)
#> b0 1.1465150 17.257565
#> b1 0.9071749  1.103594
#> b2 0.1055102  1.084782

# Should return:
#         Estimate         SE      z     p          LL         UL
# b0  1.4924887636 0.69172794  2.158 0.031  0.13672691 2.84825062
# b1  0.0005759509 0.04999884  0.012 0.991 -0.09741998 0.09857188
# b2 -1.0837844594 0.59448206 -1.823 0.068 -2.24894789 0.08137897
#     exp(Estimate)   exp(LL)   exp(UL)
# b0      4.4481522 1.1465150 17.257565
# b1      1.0005761 0.9071749  1.103594
# b2      0.3383128 0.1055102  1.084782

```
