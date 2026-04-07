# Meta-regression analysis for 2-group standardized mean differences

This function estimates the intercept and slope coefficients in a
meta-regression model where the dependent variable is a 2-group
standardized mean difference. The estimates are OLS estimates with
robust standard errors that accommodate residual heteroscedasticity. Use
the unweighted variance standardizer for 2-group experimental designs,
and use the weighted variance standardizer for 2-group nonexperimental
designs. A single-group standardizer can be used in either experimental
or nonexperimental designs.

For more details, see Section 3.4 of Bonett (2021, Volume 5).

## Usage

``` r
meta.lm.stdmean2(alpha, m1, m2, sd1, sd2, n1, n2, X, stdzr)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- m1:

  vector of estimated means for group 1

- m2:

  vector of estimated means for group 2

- sd1:

  vector of estimated SDs for group 1

- sd2:

  vector of estimated SDs for group 2

- n1:

  vector of group 1 sample sizes

- n2:

  vector of group 2 sample sizes

- X:

  matrix of predictor values

- stdzr:

  - set to 0 for square root unweighted average variance standardizer

  - set to 1 for group 1 SD standardizer

  - set to 2 for group 2 SD standardizer

  - set to 3 for square root weighted average variance standardizer

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

Bonett DG (2009). “Meta-analytic interval estimation for standardized
and unstandardized mean differences.” *Psychological Methods*,
**14**(3), 225–238. ISSN 1939-1463,
[doi:10.1037/a0016619](https://doi.org/10.1037/a0016619) .

Bonett DG (2021). *Statistical Methods for Psychologists, Vol 1-5,
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
n1 <- c(65, 30, 29, 45, 50)
n2 <- c(67, 32, 31, 20, 52)
m1 <- c(31.1, 32.3, 31.9, 29.7, 33.0)
m2 <- c(34.1, 33.2, 30.6, 28.7, 26.5)
sd1 <- c(7.1, 8.1, 7.8, 6.8, 7.6)
sd2 <- c(7.8, 7.3, 7.5, 7.2, 6.8)
x1 <- c(4, 6, 7, 7, 8)
X <- matrix(x1, 5, 1)
meta.lm.stdmean2(.05, m1, m2, sd1, sd2, n1, n2, X, 0)
#>    Estimate      SE      z p      LL      UL
#> b0  -1.6988 0.41080 -4.135 0 -2.5040 -0.8937
#> b1   0.2872 0.06498  4.419 0  0.1598  0.4145

# Should return:
#     Estimate      SE      z p      LL      UL
# b0   -1.6988 0.41080 -4.135 0 -2.5040 -0.8937
# b1    0.2872 0.06498  4.419 0  0.1598  0.4145

```
