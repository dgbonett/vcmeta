# Confidence interval for an average mean difference from paired-samples studies

Computes the estimate, standard error, and confidence interval for an
average mean difference from two or more paired-samples studies. A
Satterthwaite adjustment to the degrees of freedom is used to improve
the accuracy of the confidence interval for the average effect size.
Equality of variances within or across studies is not assumed.

For more details, see Chapter 2 of Bonett (2021, Volume 5).

## Usage

``` r
meta.ave.mean.ps(alpha, m1, m2, sd1, sd2, cor, n, bystudy = TRUE)
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

- bystudy:

  logical to also return each study estimate (TRUE) or not

## Value

Returns a matrix. The first row is the average estimate across all
studies. If bystudy is TRUE, there is 1 additional row for each study.
The matrix has the following columns:

- Estimate - estimated effect size

- SE - standard error

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

- df - degrees of freedom

## References

Bonett DG (2009). “Meta-analytic interval estimation for standardized
and unstandardized mean differences.” *Psychological Methods*,
**14**(3), 225–238. ISSN 1939-1463,
[doi:10.1037/a0016619](https://doi.org/10.1037/a0016619) .

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
meta.ave.mean.ps(.05, m1, m2, sd1, sd2, cor, n, bystudy = TRUE)
#>         Estimate        SE        LL         UL     df
#> Average    -3.25 0.2340603 -3.713965 -2.7860354 107.66
#> Study 1    -2.00 0.5672507 -3.160158 -0.8398421  29.00
#> Study 2    -2.00 0.4227434 -2.849535 -1.1504653  49.00
#> Study 3    -5.00 0.5335104 -6.091151 -3.9088487  29.00
#> Study 4    -4.00 0.3023716 -4.603215 -3.3967852  69.00

# Should return:
#         Estimate        SE        LL         UL     df
# Average    -3.25 0.2340603 -3.713965 -2.7860352 107.66
# Study 1    -2.00 0.5672507 -3.160158 -0.8398421  29.00
# Study 2    -2.00 0.4227434 -2.849535 -1.1504653  49.00
# Study 3    -5.00 0.5335104 -6.091151 -3.9088487  29.00
# Study 4    -4.00 0.3023716 -4.603215 -3.3967852  69.00

```
