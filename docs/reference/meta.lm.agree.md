# Meta-regression analysis for G agreement indices

This function estimates the intercept and slope coefficients in a
meta-regression model where the dependent variable is a G-index of
agreement. The estimates are OLS estimates with robust standard errors
that accomodate residual heteroscedasticity.

For more details, see Section 3.4 of Bonett (2021, Volume 5).

## Usage

``` r
meta.lm.agree(alpha, f11, f12, f21, f22, X)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- f11:

  vector of frequency counts in cell 1,1

- f12:

  vector of frequency counts in cell 1,2

- f21:

  vector of frequency counts in cell 2,1

- f22:

  vector of frequency counts in cell 2,2

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

Bonett DG (2022). “Statistical inference for G-indices of agreement.”
*Journal of Educational and Behavioral Statistics*, **47**(4), 438–458.
[doi:10.3102/10769986221088561](https://doi.org/10.3102/10769986221088561)
.

Bonett DG (2021). *Statistical Methods for Psychologists, Vol 1-5,
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
f11 <- c(40, 20, 25, 30)
f12 <- c(3, 2, 2, 1)
f21 <- c(7, 6, 8, 6)
f22 <- c(26, 25, 13, 25)
x1 <- c(1, 1, 4, 6)
x2 <- c(1, 1, 0, 0)
X <- matrix(cbind(x1, x2), 4, 2)
meta.lm.agree(.05, f11, f12, f21, f22, X)
#>    Estimate      SE     z     p      LL     UL
#> b0   0.1905 0.38773 0.491 0.623 -0.5695 0.9504
#> b1   0.0952 0.07142 1.334 0.182 -0.0447 0.2352
#> b2   0.4205 0.32384 1.299 0.194 -0.2142 1.0552

# Should return:
#    Estimate      SE     z     p      LL     UL
# b0   0.1905 0.38773 0.491 0.623 -0.5695 0.9504
# b1   0.0952 0.07142 1.334 0.182 -0.0447 0.2352
# b2   0.4205 0.32384 1.299 0.194 -0.2142 1.0552

```
