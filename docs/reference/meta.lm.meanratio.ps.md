# Meta-regression analysis for paired-samples log mean ratios

This function estimates the intercept and slope coefficients in a
meta-regression model where the dependent variable is a paired-samples
log mean ratio. The estimates are OLS estimates with robust standard
errors that accommodate residual heteroscedasticity. The exponentiated
slope estimate for a predictor variable describes a multiplicative
change in the mean ratio associated with a 1-unit increase in that
predictor variable, controlling for all other predictor variables in the
model.

For more details, see Section 3.4 of Bonett (2021, Volume 5).

## Usage

``` r
meta.lm.meanratio.ps(alpha, m1, m2, sd1, sd2, cor, n, X)
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

- z - z-value

- p - p-value

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

- exp(Estimate) - the exponentiated estimate

- exp(LL) - lower limit of the exponentiated confidence interval

- exp(UL) - upper limit of the exponentiated confidence interval

## References

Bonett DG, Price RM (2020). “Confidence intervals for ratios of means
and medians.” *Journal of Educational and Behavioral Statistics*,
**45**(6), 750–770. ISSN 1076-9986,
[doi:10.3102/1076998620934125](https://doi.org/10.3102/1076998620934125)
.

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
meta.lm.meanratio.ps(.05, m1, m2, sd1, sd2, cor, n, X)
#>      Estimate         SE    z     p           LL        UL exp(Estimate)
#> b0 0.50957008 0.13000068 3.92 0.000  0.254773424 0.7643667      1.664575
#> b1 0.07976238 0.04133414 1.93 0.054 -0.001251047 0.1607758      1.083030
#>      exp(LL)  exp(UL)
#> b0 1.2901693 2.147634
#> b1 0.9987497 1.174422

# Should return: 
#      Estimate         SE           LL        UL     z     p
# b0 0.50957008 0.13000068  0.254773424 0.7643667 3.920 0.000
# b1 0.07976238 0.04133414 -0.001251047 0.1607758 1.920 0.054
#     exp(Estimate)   exp(LL)  exp(UL)
# b0       1.664575 1.2901693 2.147634
# b1       1.083030 0.9987497 1.174422

```
