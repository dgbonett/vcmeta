# Confidence interval for a linear contrast of standardized mean differences from 2-group studies

Computes the estimate, standard error, and confidence interval for a
linear contrast of 2-group standardized mean differences from two or
more studies. Equality of variances within or across studies is not
assumed. Use the square root average variance standardizer (stdzr = 0)
for 2-group experimental designs. Use the square root weighted variance
standardizer (stdzr = 3) for 2-group nonexperimental designs with simple
random sampling. The stdzr = 1 and stdzr = 2 options can be used with
either experimental or nonexperimental designs.

For more details, see Section 3.2 of Bonett (2021, Volume 5).

## Usage

``` r
meta.lc.stdmean2(alpha, m1, m2, sd1, sd2, n1, n2, v, stdzr)
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

- v:

  vector of contrast coefficients

- stdzr:

  - set to 0 for square root unweighted average variance standardizer

  - set to 1 for group 1 SD standardizer

  - set to 2 for group 2 SD standardizer

  - set to 3 for square root weighted average variance standardizer

## Value

Returns 1-row matrix with the following columns:

- Estimate - estimated linear contrast

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
m1 <- c(45.1, 39.2, 36.3, 34.5)
m2 <- c(30.0, 35.1, 35.3, 36.2)
sd1 <- c(10.7, 10.5, 9.4, 11.5)
sd2 <- c(12.3, 12.0, 10.4, 9.6)
n1 <- c(40, 20, 50, 25)
n2 <- c(40, 20, 48, 26)
v <- c(.5, .5, -.5, -.5)
meta.lc.stdmean2(.05, m1, m2, sd1, sd2, n1, n2, v, 0)
#>          Estimate      SE     LL     UL
#> Contrast   0.8558 0.27092 0.3248 1.3868

# Should return: 
#          Estimate      SE     LL     UL
# Contrast   0.8558 0.27092 0.3248 1.3868

```
