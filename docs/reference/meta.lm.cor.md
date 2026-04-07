# Meta-regression analysis for Pearson or partial correlations

This function estimates the intercept and slope coefficients in a
meta-regression model where the dependent variable is a
Fisher-transformed Pearson or partial correlation. The estimates are OLS
estimates with robust standard errors that accommodate residual
heteroscedasticity. The correlations are Fisher-transformed and hence
the parameter estimates do not have a simple interpretation. However,
the hypothesis test results can be used to decide if a population slope
is either positive or negative.

For more details, see Section 3.4 of Bonett (2021, Volume 5).

## Usage

``` r
meta.lm.cor(alpha, n, cor, s, X)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- n:

  vector of sample sizes

- cor:

  vector of estimated Pearson or partial correlations

- s:

  number of control variables

- X:

  matrix of predictor values

## Value

Returns a matrix. The first row is for the intercept with one additional
row per predictor. The matrix has the following columns:

- Estimate - OLS estimate

- SE - Standard error

- z - z-value

- p - p-value

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## References

Bonett DG (2021). *Statistical Methods for Psychologists, Vol 1-5,
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
n <- c(55, 190, 65, 35)
cor <- c(.40, .65, .60, .45)
q <- 0
x1 <- c(18, 25, 23, 19)
X <- matrix(x1, 4, 1)
meta.lm.cor(.05, n, cor, q, X)
#>    Estimate      SE      z     p      LL     UL
#> b0  -0.4783 0.48632 -0.984 0.325 -1.4315 0.4748
#> b1   0.0505 0.02128  2.371 0.018  0.0088 0.0922

# Should return: 
#    Estimate      SE      z     p      LL     UL
# b0  -0.4783 0.48632 -0.984 0.325 -1.4315 0.4748
# b1   0.0505 0.02128  2.371 0.018  0.0088 0.0922

```
