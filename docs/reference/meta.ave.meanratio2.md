# Confidence interval for an average mean ratio from 2-group studies

Computes the estimate, standard error, and confidence interval for a
geometric average mean ratio from two or more 2-group studies. A
Satterthwaite adjustment to the degrees of freedom is used to improve
the accuracy of the confidence intervals. Equality of variances within
or across studies is not assumed.

For more details, see Chapter 2 of Bonett (2021, Volume 5).

## Usage

``` r
meta.ave.meanratio2(alpha, m1, m2, sd1, sd2, n1, n2, bystudy = TRUE)
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

- exp(Estimate) - exponentiated estimate

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
m1 <- c(7.4, 6.9)
m2 <- c(6.3, 5.7)
sd1 <- c(1.7, 1.5)
sd2 <- c(2.3, 2.0)
n1 <- c(40, 20)
n2 <- c(40, 20)
meta.ave.meanratio2(.05, m1, m2, sd1, sd2, n1, n2, bystudy = TRUE)
#>          Estimate         SE          LL        UL exp(Estimate)  exp(LL)
#> Average 0.1759928 0.05738065 0.061437026 0.2905486      1.192429 1.063364
#> Study 1 0.1609304 0.06820167 0.024749712 0.2971110      1.174603 1.025059
#> Study 2 0.1910552 0.09229675 0.002986265 0.3791242      1.210526 1.002991
#>          exp(UL)       df
#> Average 1.337161 66.26000
#> Study 1 1.345965 65.69929
#> Study 2 1.461004 31.71341

# Should return:
#          Estimate         SE          LL        UL  exp(Estimate)
# Average 0.1759928 0.05738065 0.061437186 0.2905484       1.192429
# Study 1 0.1609304 0.06820167 0.024749712 0.2971110       1.174603
# Study 2 0.1910552 0.09229675 0.002986265 0.3791242       1.210526
#          exp(LL)  exp(UL)    df
# Average 1.063364 1.337161 66.27
# Study 1 1.025059 1.345965 65.70
# Study 2 1.002991 1.461004 31.71

```
