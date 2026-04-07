# Meta-regression analysis for semipartial correlations

This function estimates the intercept and slope coefficients in a
meta-regression model where the dependent variable is a
Fisher-transformed semipartial correlation. The estimates are OLS
estimates with robust standard errors that accommodate residual
heteroscedasticity. The correlations are Fisher-transformed and hence
the parameter estimates do not have a simple interpretation. However,
the hypothesis test results can be used to decide if a population slope
is either positive or negative.

For more details, see Section 3.4 of Bonett (2021, Volume 5).

## Usage

``` r
meta.lm.semipart(alpha, n, cor, r2, X)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- n:

  vector of sample sizes

- cor:

  vector of estimated semipartial correlations

- r2:

  vector of squared multiple correlations for a model that includes the
  IV and all control variables

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
n <- c(128, 97, 210, 217)
cor <- c(.35, .41, .44, .39)
r2 <- c(.29, .33, .36, .39)
x1 <- c(18, 25, 23, 19)
X <- matrix(x1, 4, 1)
meta.lm.semipart(.05, n, cor, r2, X)
#>    Estimate      SE     z     p      LL     UL
#> b0   0.1970 0.30618 0.643 0.520 -0.4031 0.7971
#> b1   0.0106 0.01457 0.725 0.468 -0.0180 0.0391

# Should return: 
#    Estimate      SE     z     p      LL     UL
# b0   0.1970 0.30618 0.643 0.520 -0.4031 0.7971
# b1   0.0106 0.01457 0.725 0.468 -0.0180 0.0391

```
