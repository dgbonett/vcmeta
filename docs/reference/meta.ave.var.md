# Confidence interval for an average variance

Computes the estimate and confidence interval for an average variance
from two or more studies. The estimated average variance or the upper
confidence limit could be used as a variance planning value in sample
size planning. The confidence intervals assume normality, and this
function is not recommended if the variances have been estimated from
leptokurtic data.

For more details, see Chapter 2 of Bonett (2021, Volume 5).

## Usage

``` r
meta.ave.var(alpha, var, n, bystudy = TRUE)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- var:

  vector of sample variances

- n:

  vector of sample sizes

- bystudy:

  logical to also return each study estimate (TRUE) or not

## Value

Returns a matrix. The first row is the average estimate across all
studies. If bystudy is TRUE, there is 1 additional row for each study.
The matrix has the following columns:

- Estimate - estimated variance

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## References

Bonett DG (2021). *Statistical Methods for Psychologists, Vol 1-5,
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
var <- c(26.63, 22.45, 34.12)
n <- c(40, 30, 50)
meta.ave.var(.05, var, n, bystudy = TRUE)
#>         Estimate       LL       UL
#> Average 27.73333 21.45679 35.84589
#> Study 1 26.63000 17.86939 43.90614
#> Study 2 22.45000 14.23923 40.57127
#> Study 3 34.12000 23.80835 52.98319

# Should return:
#         Estimate       LL       UL
# Average 27.73333 21.45679 35.84589
# Study 1 26.63000 17.86939 43.90614
# Study 2 22.45000 14.23923 40.57127
# Study 3 34.12000 23.80835 52.98319

```
