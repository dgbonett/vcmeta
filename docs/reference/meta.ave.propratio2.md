# Confidence interval for an average proportion ratio from 2-group studies

Computes the estimate, standard error, and confidence interval for a
geometric average proportion ratio from two or more studies.

For more details, see Chapter 2 of Bonett (2021, Volume 5).

## Usage

``` r
meta.ave.propratio2(alpha, f1, f2, n1, n2, bystudy = TRUE)
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

## References

Price RM, Bonett DG (2008). “Confidence intervals for a ratio of two
independent binomial proportions.” *Statistics in Medicine*, **27**(26),
5497–5508. ISSN 02776715,
[doi:10.1002/sim.3376](https://doi.org/10.1002/sim.3376) .

Bonett DG (2021). *Statistical Methods for Psychologists, Vol 1-5,
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
n1 <- c(204, 201, 932, 130, 77)
n2 <- c(106, 103, 415, 132, 83)
f1 <- c(24, 40, 93, 14, 5)
f2 <- c(12, 9, 28, 3, 1)
meta.ave.propratio2(.05, f1, f2, n1, n2, bystudy = TRUE)
#>           Estimate        SE          LL        UL exp(Estimate)   exp(LL)
#> Average 0.84705608 0.2528742  0.35143178 1.3426804      2.332769 1.4211008
#> Study 1 0.03604257 0.3297404 -0.61023681 0.6823220      1.036700 0.5432222
#> Study 2 0.81008932 0.3442007  0.13546839 1.4847103      2.248109 1.1450730
#> Study 3 0.38746839 0.2065227 -0.01730864 0.7922454      1.473246 0.9828403
#> Study 4 1.49316811 0.6023296  0.31262374 2.6737125      4.451175 1.3670071
#> Study 5 1.50851199 0.9828420 -0.41782290 3.4348469      4.520000 0.6584788
#>           exp(UL)
#> Average  3.829294
#> Study 1  1.978466
#> Study 2  4.413686
#> Study 3  2.208350
#> Study 4 14.493677
#> Study 5 31.026662

#  Should return:
#           Estimate        SE          LL        UL 
# Average 0.84705608 0.2528742  0.35143178 1.3426804
# Study 1 0.03604257 0.3297404 -0.61023681 0.6823220
# Study 2 0.81008932 0.3442007  0.13546839 1.4847103
# Study 3 0.38746839 0.2065227 -0.01730864 0.7922454
# Study 4 1.49316811 0.6023296  0.31262374 2.6737125
# Study 5 1.50851199 0.9828420 -0.41782290 3.4348469
#     exp(Estimate)   exp(LL)   exp(UL)
# Average  2.332769 1.4211008  3.829294
# Study 1  1.036700 0.5432222  1.978466
# Study 2  2.248109 1.1450730  4.413686
# Study 3  1.473246 0.9828403  2.208350
# Study 4  4.451175 1.3670071 14.493677
# Study 5  4.520000 0.6584788 31.026662

```
