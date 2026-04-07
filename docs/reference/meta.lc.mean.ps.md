# Confidence interval for a linear contrast of mean differences from paired-samples studies

Computes the estimate, standard error, and confidence interval for a
linear contrast of paired-samples mean differences from two or more
studies. A Satterthwaite adjustment to the degrees of freedom is used to
improve the accuracy of the confidence interval. Equality of variances
within or across studies is not assumed.

For more details, see Section 3.2 of Bonett (2021, Volume 5).

## Usage

``` r
meta.lc.mean.ps(alpha, m1, m2, sd1, sd2, cor, n, v)
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

- Estimate - estimated linear contrast

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
v <- c(.5, .5, -.5, -.5)
meta.lc.mean.ps(.05, m1, m2, sd1, sd2, cor, n, v)
#>          Estimate        SE       LL       UL     df
#> Contrast      2.5 0.4681205 1.572071 3.427929 107.66

# Should return:
#          Estimate        SE      LL      UL     df
# Contrast      2.5 0.4681205 1.57207 3.42793 107.66

```
