# Meta-regression analysis for 2-group log mean ratios

This function estimates the intercept and slope coefficients in a
meta-regression model where the dependent variable is a 2-group log mean
ratio. The estimates are OLS estimates with robust standard errors that
accommodate residual heteroscedasticity. The exponentiated slope
estimate for a predictor variable describes a multiplicative change in
the mean ratio associated with a 1-unit increase in that predictor
variable, controlling for all other predictor variables in the model.

For more details, see Section 3.4 of Bonett (2021, Volume 5).

## Usage

``` r
meta.lm.meanratio2(alpha, m1, m2, sd1, sd2, n1, n2, X)
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

Bonett DG, Price RM (2020). “Confidence intervals for ratios of means
and medians.” *Journal of Educational and Behavioral Statistics*,
**45**(6), 750–770. ISSN 1076-9986,
[doi:10.3102/1076998620934125](https://doi.org/10.3102/1076998620934125)
.

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
meta.lm.meanratio2(.05, m1, m2, sd1, sd2, n1, n2, X)
#>       Estimate         SE      z p          LL          UL exp(Estimate)
#> b0 -0.40208954 0.09321976 -4.313 0 -0.58479692 -0.21938216     0.6689208
#> b1  0.06831545 0.01484125  4.603 0  0.03922712  0.09740377     1.0707030
#>     exp(LL)   exp(UL)
#> b0 0.557219 0.8030148
#> b1 1.040007 1.1023054

# Should return:
#       Estimate         SE          LL          UL      z p
# b0 -0.40208954 0.09321976 -0.58479692 -0.21938216 -4.313 0
# b1  0.06831545 0.01484125  0.03922712  0.09740377  4.603 0
#    exp(Estimate)  exp(LL)   exp(UL)
# b0     0.6689208 0.557219 0.8030148
# b1     1.0707030 1.040007 1.1023054

```
