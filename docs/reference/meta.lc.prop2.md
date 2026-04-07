# Confidence interval for a linear contrast of proportion differences in 2-group studies

Computes the estimate, standard error, and adjusted Wald confidence
interval for a linear contrast of 2-group proportion differences from
two or more studies.

For more details, see Section 3.2 of Bonett (2021, Volume 5).

## Usage

``` r
meta.lc.prop2(alpha, f1, f2, n1, n2, v)
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

- v:

  vector of contrast coefficients

## Value

Returns 1-row matrix with the following columns:

- Estimate - estimated linear contrast

- SE - standard error

- LL - lower limit of the adjusted Wald confidence interval

- UL - upper limit of the adjusted Wald confidence interval

## References

Bonett DG, Price RM (2014). “Meta-analysis methods for risk
differences.” *British Journal of Mathematical and Statistical
Psychology*, **67**(3), 371–387. ISSN 00071102,
[doi:10.1111/bmsp.12024](https://doi.org/10.1111/bmsp.12024) .

Bonett DG (2021). *Statistical Methods for Psychologists, Vol 1-5,
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
n1 <- c(50, 150, 150)
n2 <- c(50, 150, 150)
f1 <- c(16, 50, 25)
f2 <- c(7, 15, 20)
v <- c(1, -1, 0)
meta.lc.prop2(.05, f1, f2, n1, n2, v)
#>             Estimate         SE         LL        UL
#> Contrast -0.05466931 0.09401019 -0.2389259 0.1295873

# Should return:
#             Estimate         SE         LL        UL
# Contrast -0.05466931 0.09401019 -0.2389259 0.1295873

```
