# Confidence interval for a linear contrast of G-index coefficients

Computes the estimate, standard error, and adjusted Wald confidence
interval for a linear contrast of G-index of agreement coefficients from
two or more studies. This function assumes that two raters each provide
a dichotomous rating for a sample of objects.

For more details, see Setion 3.2 of Bonett (2021, Volume 5).

## Usage

``` r
meta.lc.agree(alpha, f11, f12, f21, f22, v)
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

Bonett DG (2022). “Statistical inference for G-indices of agreement.”
*Journal of Educational and Behavioral Statistics*, **47**(4), 438–458.
[doi:10.3102/10769986221088561](https://doi.org/10.3102/10769986221088561)
.

Bonett DG (2021). *Statistical Methods for Psychologists, Vol 1-5,
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
f11 <- c(43, 56, 49)
f12 <- c(7, 2, 9)
f21 <- c(3, 5, 5)
f22 <- c(37, 54, 39)
v <- c(.5, .5, -1)
meta.lc.agree(.05, f11, f12, f21, f22, v)
#>          Estimate      SE     LL     UL
#> Contrast   0.1023 0.07972 -0.054 0.2585

# Should return:
#          Estimate      SE      L     UL
# Contrast   0.1023 0.07972 -0.054 0.2585

```
