# Meta-regression analysis for correlations

This function estimates the intercept and slope coefficients in a
meta-regression model where the dependent variable is a
Fisher-transformed correlation. The correlations can be of different
types (e.g., Pearson, partial, Spearman). The estimates are OLS
estimates with robust standard errors that accommodate residual
heteroscedasticity. This function uses estimated correlations and their
standard errors as input. The correlations are Fisher-transformed and
hence the parameter estimates do not have a simple interpretation.
However, the hypothesis test results can be used to decide if a
population slope is either positive or negative.

For more details, see Section 3.4 of Bonett (2021, Volume 5).

## Usage

``` r
meta.lm.cor.gen(alpha, cor, se, X)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- cor:

  vector of estimated correlations

- se:

  number of control variables

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
cor <- c(.40, .65, .60, .45)
se <- c(.182, .114, .098, .132)
x1 <- c(18, 25, 23, 19)
X <- matrix(x1, 4, 1)
meta.lm.cor.gen(.05, cor, se, X)
#>    Estimate      SE      z     p
#> b0  -0.4783 0.63428 -0.754 0.451
#> b1   0.0505 0.02880  1.753 0.080

# Should return: 
#    Estimate      SE      z     p
# b0  -0.4783 0.63428 -0.754 0.451
# b1   0.0505 0.02880  1.753 0.080

```
