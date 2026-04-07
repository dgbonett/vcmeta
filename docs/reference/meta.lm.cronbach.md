# Meta-regression analysis for Cronbach reliabilities

This function estimates the intercept and slope coefficients in a
meta-regression model where the dependent variable is a log-complement
Cronbach reliablity. The estimates are OLS estimates with robust
standard errors that accommodate residual heteroscedasticity. The
exponentiated slope estimate for a predictor variable describes a
multiplicative change in non-reliability associated with a 1-unit
increase in that predictor variable, controlling for all other predictor
variables in the model.

For more details, see Section 3.4 of Bonett (2021, Volume 5).

## Usage

``` r
meta.lm.cronbach(alpha, n, rel, r, X)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- n:

  vector of sample sizes

- rel:

  vector of estimated reliabilities

- r:

  number of measurements (e.g., items)

- X:

  matrix of predictor values

## Value

Returns a matrix. The first row is for the intercept with one additional
row per predictor. The matrix has the following columns:

- Estimate - exponentiated OLS estimate

- SE - standard error

- z - z-value

- p - p-value

- LL - lower limit of the exponentiated confidence interval

- UL - upper limit of the exponentiated confidence interval

## References

Bonett DG (2010). “Varying coefficient meta-analytic methods for alpha
reliability.” *Psychological Methods*, **15**(4), 368–385. ISSN
1939-1463, [doi:10.1037/a0020142](https://doi.org/10.1037/a0020142) .

Bonett DG, Wright TA (2015). “Cronbach's alpha reliability: Interval
estimation, hypothesis testing, and sample size planning.” *Journal of
Organizational Behavior*, **36**(1), 3–15. ISSN 08943796,
[doi:10.1002/job.1960](https://doi.org/10.1002/job.1960) .

Bonett DG (2021). *Statistical Methods for Psychologists, Vol 1-5,
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
n <- c(583, 470, 546, 680)
rel <- c(.91, .89, .90, .89)
x1 <- c(1, 0, 0, 0)
X <- matrix(x1, 4, 1)
meta.lm.cronbach(.05, n, rel, 10, X)
#>    Estimate      SE       z     p      LL      UL
#> b0  -2.2408 0.03676 -60.960 0.000 -2.3129 -2.1688
#> b1  -0.1689 0.07205  -2.344 0.019 -0.3101 -0.0277

# Should return:
#    Estimate      SE       z     p      LL      UL
# b0  -2.2408 0.03676 -60.960 0.000 -2.3129 -2.1688
# b1  -0.1689 0.07205  -2.344 0.019 -0.3101 -0.0277

```
