# Meta-regression analysis for odds ratios

This function estimates the intercept and slope coefficients in a
meta-regression model where the dependent variable is a log odds ratio.
The estimates are OLS estimates with robust standard errors that
accommodate residual heteroscedasticity. The exponentiated slope
estimate for a predictor variable describes a multiplicative change in
the odds ratio associated with a 1-unit increase in that predictor
variable, controlling for all other predictor variables in the model.

For more details, see Section 3.4 of Bonett (2021, Volume 5).

## Usage

``` r
meta.lm.oddsratio(alpha, f1, f2, n1, n2, X)
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

Bonett DG, Price RM (2015). “Varying coefficient meta-analysis methods
for odds ratios and risk ratios.” *Psychological Methods*, **20**(3),
394–406. ISSN 1939-1463,
[doi:10.1037/met0000032](https://doi.org/10.1037/met0000032) .

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
meta.lm.oddsratio(.05, f1, f2, n1, n2, X)
#>        Estimate         SE      z     p         LL         UL exp(Estimate)
#> b0  1.541895013 0.69815801  2.209 0.027  0.1735305 2.91025958     4.6734381
#> b1 -0.004417932 0.04840623 -0.091 0.927 -0.0992924 0.09045653     0.9955918
#> b2 -1.071122269 0.60582695 -1.768 0.077 -2.2585213 0.11627674     0.3426238
#>      exp(LL)   exp(UL)
#> b0 1.1894969 18.361564
#> b1 0.9054779  1.094674
#> b2 0.1045049  1.123307

# Should return:
#        Estimate         SE      z     p         LL         UL
# b0  1.541895013 0.69815801  2.209 0.027  0.1735305 2.91025958
# b1 -0.004417932 0.04840623 -0.091 0.927 -0.0992924 0.09045653
# b2 -1.071122269 0.60582695 -1.768 0.077 -2.2585213 0.11627674
#    exp(Estimate)   exp(LL)   exp(UL)
# b0     4.6734381 1.1894969 18.361564
# b1     0.9955918 0.9054779  1.094674
# b2     0.3426238 0.1045049  1.123307

```
