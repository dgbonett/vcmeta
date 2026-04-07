# Confidence interval for an average mean difference from 2-group studies

Computes the estimate, standard error, and confidence interval for an
average mean difference from two or more 2-group studies. A
Satterthwaite adjustment to the degrees of freedom is used to improve
the accuracy of the confidence intervals. Equality of variances within
or across studies is not assumed.

For more details, see Chapter 2 of Bonett (2021, Volume 5).

## Usage

``` r
meta.ave.mean2(alpha, m1, m2, sd1, sd2, n1, n2, bystudy = TRUE)
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
m1 <- c(7.4, 6.9)
m2 <- c(6.3, 5.7)
sd1 <- c(1.72, 1.53)
sd2 <- c(2.35, 2.04)
n1 <- c(40, 60)
n2 <- c(40, 60)
meta.ave.mean2(.05, m1, m2, sd1, sd2, n1, n2, bystudy = TRUE)
#>         Estimate        SE        LL       UL        df
#> Average     1.15 0.2830183 0.5904369 1.709563 139.41000
#> Study 1     1.10 0.4604590 0.1819748 2.018025  71.46729
#> Study 2     1.20 0.3292036 0.5475574 1.852443 109.42136

# Should return:
#         Estimate        SE        LL       UL     df
# Average     1.15 0.2830183 0.5904369 1.709563 139.41
# Study 1     1.10 0.4604590 0.1819748 2.018025  71.47
# Study 2     1.20 0.3292036 0.5475574 1.852443 109.42

```
