# Meta-regression analysis for paired-samples proportion differences

This function estimates the intercept and slope coefficients in a
meta-regression model where the dependent variable is a paired-samples
proportion difference. The estimates are OLS estimates with robust
standard errors that accommodate residual heteroscedasticity.

For more details, see Section 3.4 of Bonett (2021, Volume 5).

## Usage

``` r
meta.lm.prop.ps(alpha, f11, f12, f21, f22, X)
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

Bonett DG, Price RM (2012). “Adjusted wald confidence interval for a
difference of binomial proportions based on paired data.” *Journal of
Educational and Behavioral Statistics*, **37**(4), 479–488. ISSN
1076-9986,
[doi:10.3102/1076998611411915](https://doi.org/10.3102/1076998611411915)
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
meta.lm.prop.ps(.05, f11, f12, f21, f22, X)
#>       Estimate         SE      z     p          LL         UL
#> b0 -0.21113402 0.21119823 -1.000 0.317 -0.62507494 0.20280690
#> b1  0.02185567 0.03861947  0.566 0.571 -0.05383711 0.09754845
#> b2  0.12575138 0.17655623  0.712 0.476 -0.22029248 0.47179524

# Should return: 
#       Estimate         SE      z     p          LL         UL
# b0 -0.21113402 0.21119823 -1.000 0.317 -0.62507494 0.20280690
# b1  0.02185567 0.03861947  0.566 0.571 -0.05383711 0.09754845
# b2  0.12575138 0.17655623  0.712 0.476 -0.22029248 0.47179524

```
