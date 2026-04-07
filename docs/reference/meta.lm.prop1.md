# Meta-regression analysis for 1-group proportions

This function estimates the intercept and slope coefficients in a
meta-regression model where the dependent variable is a proportion from
one group. The estimates are OLS estimates with robust standard errors
that accomodate residual heteroscedasticity.

For more details, see Section 3.4 of Bonett (2021, Volume 5).

## Usage

``` r
meta.lm.prop1(alpha, f, n, X)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- f:

  vector of frequency counts

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

## References

Bonett DG (2021). *Statistical Methods for Psychologists, Vol 1-5,
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
f <- c(38, 26, 24, 15, 45, 38)
n <- c(80, 60, 70, 50, 180, 200)
x1 <- c(10, 15, 18, 22, 24, 30)
X <- matrix(x1, 6, 1)
meta.lm.prop1(.05, f, n, X)
#>       Estimate         SE      z p          LL           UL
#> b0  0.63262816 0.06845707  9.241 0  0.49845477  0.766801546
#> b1 -0.01510565 0.00290210 -5.205 0 -0.02079367 -0.009417641

# Should return: 
#       Estimate         SE      z p          LL           UL
# b0  0.63262816 0.06845707  9.241 0  0.49845477  0.766801546
# b1 -0.01510565 0.00290210 -5.205 0 -0.02079367 -0.009417641

```
