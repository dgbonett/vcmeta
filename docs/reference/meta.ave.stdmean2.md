# Confidence interval for an average standardized mean difference from 2-group studies

Computes the estimate, standard error, and confidence interval for an
average standardized mean difference from two or more 2-group studies.
Square root unweighted variances, square root weighted variances, and
single group standard deviation are options for the standardizer.
Equality of variances within or across studies is not assumed.

For more details, see Chapter 2 of Bonett (2021, Volume 5).

## Usage

``` r
meta.ave.stdmean2(alpha, m1, m2, sd1, sd2, n1, n2, stdzr, bystudy = TRUE)
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

- stdzr:

  - set to 0 for square root unweighted average variance standardizer

  - set to 1 for group 1 SD standardizer

  - set to 2 for group 2 SD standardizer

  - set to 3 for square root weighted average variance standardizer

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

## References

Bonett DG (2009). “Meta-analytic interval estimation for standardized
and unstandardized mean differences.” *Psychological Methods*,
**14**(3), 225–238. ISSN 1939-1463,
[doi:10.1037/a0016619](https://doi.org/10.1037/a0016619) .

Bonett DG (2021). *Statistical Methods for Psychologists, Vol 1-5,
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
m1 <- c(21.9, 23.1, 19.8)
m2 <- c(16.1, 17.4, 15.0)
sd1 <- c(3.82, 3.95, 3.67)
sd2 <- c(3.21, 3.30, 3.02)
n1 <- c(40, 30, 24)
n2 <- c(40, 28, 25)
meta.ave.stdmean2(.05, m1, m2, sd1, sd2, n1, n2, 0, bystudy = TRUE)
#>         Estimate      SE     LL     UL
#> Average   1.5261 0.17343 1.1862 1.8661
#> Study 1   1.6439 0.26290 1.1286 2.1592
#> Study 2   1.5661 0.30563 0.9671 2.1652
#> Study 3   1.4283 0.32892 0.7836 2.0729

# Should return: 
#         Estimate      SE     LL     UL
# Average   1.5261 0.17343 1.1862 1.8661
# Study 1   1.6439 0.26290 1.1286 2.1592
# Study 2   1.5661 0.30563 0.9671 2.1652
# Study 3   1.4283 0.32892 0.7836 2.0729

m1 <- c(41.2, 43.2, 49.1, 40.8)
m2 <- c(36.4, 37.1, 35.9, 31.4)
sd1 <- c(4.92, 4.75, 4.87, 5.01)
sd2 <- c(4.35, 4.24, 4.12, 4.87)
n1 <- c(42, 58, 62, 39)
n2 <- c(67, 70, 84, 45)
meta.ave.stdmean2(.05, m1, m2, sd1, sd2, n1, n2, 3, bystudy = FALSE)
#>         Estimate      SE     LL     UL
#> Average   1.8078 0.11648 1.5795 2.0361

# Should return: 
#         Estimate      SE     LL     UL
# Average   1.8078 0.11648 1.5795 2.0361

```
