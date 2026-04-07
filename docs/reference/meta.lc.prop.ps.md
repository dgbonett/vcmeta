# Confidence interval for a linear contrast of proportion differences in paired-samples studies

Computes the estimate, standard error, and adjusted Wald confidence
interval for a linear contrast of paired-samples proportion differences
from two or more studies.

For more details, see Setion 3.2 of Bonett (2021, Volume 5).

## Usage

``` r
meta.lc.prop.ps(alpha, f11, f12, f21, f22, v)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- f11:

  vector of frequency counts in cell 1,1

- f12:

  vector of frequency counts in cell 1,2

- f21:

  vector of frequency counts in cell 2,1

- f22:

  vector of frequency counts in cell 2,2

- v:

  vector of contrast coefficients

## Value

Returns 1-row matrix with the following columns:

- Estimate - estimated linear contrast

- SE - standard error

- LL - lower limit of the adjusted Wald confidence interval

- UL - upper limit of the adjusted Wald confidence interval

## References

Bonett DG, Price RM (2012). “Adjusted wald confidence interval for a
difference of binomial proportions based on paired data.” *Journal of
Educational and Behavioral Statistics*, **37**(4), 479–488. ISSN
1076-9986,
[doi:10.3102/1076998611411915](https://doi.org/10.3102/1076998611411915)
.

Bonett DG (2021). *Statistical Methods for Psychologists, Vol 1-5,
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
f11 <- c(17, 28, 19)
f12 <- c(43, 56, 49)
f21 <- c(3, 5, 5)
f22 <- c(37, 54, 39)
v <- c(.5, .5, -1)
meta.lc.prop.ps(.05, f11, f12, f21, f22, v)
#>             Estimate         SE         LL       UL
#> Contrast -0.01436285 0.06511285 -0.1419817 0.113256

# Should return:
#              Estimate         SE         LL       UL
#  Contrast -0.01436285 0.06511285 -0.1419817 0.113256

```
