# Confidence interval for a log-linear contrast of mean ratios from 2-group studies

Computes the estimate, standard error, and confidence interval for a
log-linear contrast of 2-group mean ratios from two or more studies. A
Satterthwaite adjustment to the degrees of freedom is used to improve
the accuracy of the confidence interval. Equality of variances within or
across studies is not assumed.

For more details, see Section 3.2 of Bonett (2021, Volume 5).

## Usage

``` r
meta.lc.meanratio2(alpha, m1, m2, sd1, sd2, n1, n2, v)
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
m1 <- c(45.1, 39.2, 36.3, 34.5)
m2 <- c(30.0, 35.1, 35.3, 36.2)
sd1 <- c(10.7, 10.5, 9.4, 11.5)
sd2 <- c(12.3, 12.0, 10.4, 9.6)
n1 <- c(40, 20, 50, 25)
n2 <- c(40, 20, 48, 26)
v <- c(.5, .5, -.5, -.5)
meta.lc.meanratio2(.05, m1, m2, sd1, sd2, n1, n2, v)
#>           Estimate         SE        LL        UL exp(Estimate)  exp(LL)
#> Contrast 0.2691627 0.07959269 0.1119191 0.4264064      1.308868 1.118422
#>           exp(UL)     df
#> Contrast 1.531743 152.87

# Should return:
#           Estimate         SE        LL        UL  exp(Estimate)
# Contrast 0.2691627 0.07959269 0.1119191 0.4264064      1.308868
#           exp(LL)  exp(UL)     df
# Contrast 1.118422 1.531743 152.87

```
