# Meta-regression analysis for Spearman correlations

This function estimates the intercept and slope coefficients in a
meta-regression model where the dependent variable is a
Fisher-transformed Spearman correlation. The estimates are OLS estimates
with robust standard errors that accommodate residual
heteroscedasticity. The correlations are Fisher-transformed and hence
the parameter estimates do not have a simple interpretation. However,
the hypothesis test results can be used to decide if a population slope
is either positive or negative.

For more details, see Section 3.4 of Bonett (2021, Volume 5).

## Usage

``` r
meta.lm.spear(alpha, n, cor, X)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- n:

  vector of sample sizes

- cor:

  vector of estimated Spearman correlations

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
n <- c(150, 200, 300, 200, 350)
cor <- c(.14, .29, .16, .21, .23)
x1 <- c(18, 25, 23, 19, 24)
X <- matrix(x1, 5, 1)
meta.lm.spear(.05, n, cor, X)
#>    Estimate      SE      z     p      LL     UL
#> b0  -0.0892 0.26686 -0.334 0.738 -0.6122 0.4338
#> b1   0.0137 0.01190  1.152 0.249 -0.0096 0.0370

# Should return: 
#    Estimate      SE      z     p      LL     UL
# b0  -0.0892 0.26686 -0.334 0.738 -0.6122 0.4338
# b1   0.0137 0.01190  1.152 0.249 -0.0096 0.0370

```
