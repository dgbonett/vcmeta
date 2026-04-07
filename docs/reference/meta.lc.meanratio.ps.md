# Confidence interval for a log-linear contrast of mean ratios from paired-samples studies

Computes the estimate, standard error, and confidence interval for a
log-linear contrast of paired-sample mean ratios from two or more
studies. A Satterthwaite adjustment to the degrees of freedom is used to
improve the accuracy of the confidence interval. Equality of variances
within or across studies is not assumed.

For more details, see Section 3.2 of Bonett (2021, Volume 5).

## Usage

``` r
meta.lc.meanratio.ps(alpha, m1, m2, sd1, sd2, cor, n, v)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- m1:

  vector of estimated means for measurement 1

- m2:

  vector of estimated means for measurement 2

- sd1:

  vector of estimated SDs for measurement 1

- sd2:

  vector of estimated SDs for measurement 2

- cor:

  vector of estimated correlations for paired measurements

- n:

  vector of sample sizes

- v:

  vector of contrast coefficients

## Value

Returns 1-row matrix with the following columns:

- Estimate - estimated log-linear contrast

- SE - standard error of log-linear contrast

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

- exp(Estimate) - exponentiated log-linear contrast

- exp(LL) - lower limit of the exponentiated confidence interval

- exp(UL) - upper limit of the exponentiated confidence interval

- df - degrees of freedom

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
m1 <- c(53, 60, 53, 57)
m2 <- c(55, 62, 58, 61)
sd1 <- c(4.1, 4.2, 4.5, 4.0)
sd2 <- c(4.2, 4.7, 4.9, 4.8)
cor <- c(.72, .78, .81, .85)
n <- c(30, 50, 30, 70)
v <- c(.5, .5, -.5, -.5)
meta.lc.meanratio.ps(.05, m1, m2, sd1, sd2, cor, n, v)
#>           Estimate       SE         LL         UL exp(Estimate)  exp(LL)
#> Contrast 0.0440713 0.008265 0.02767047 0.06047213      1.045057 1.028057
#>           exp(UL)    df
#> Contrast 1.062338 98.38

# Should return:
#           Estimate       SE         LL         UL exp(Estimate)
# Contrast 0.0440713 0.008265 0.02767047 0.06047213      1.045057
#           exp(LL)  exp(UL)    df
# Contrast 1.028057 1.062338 98.38

```
