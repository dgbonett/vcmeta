# Confidence interval for a linear contrast of proportions

Computes the estimate, standard error, and an adjusted Wald confidence
interval for a linear contrast of proportions from two or more studies.

For more details, see Setion 3.2 of Bonett (2021, Volume 5).

## Usage

``` r
meta.lc.prop1(alpha, f, n, v)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- f:

  vector of frequency counts

- n:

  vector of sample sizes

- v:

  vector of contrast coefficients

## Value

Returns 1-row matrix with the following columns:

- Estimate - estimated linear contrast

- SE - standard error

- LL - lower limit of the adjusted Wald confidence interval

- UL - upper limit of the adjusted Wald confidence interval

## References

Price RM, Bonett DG (2004). “An improved confidence interval for a
linear function of binomial proportions.” *Computational Statistics and
Data Analysis*, **45**(3), 449–456. ISSN 01679473,
[doi:10.1016/S0167-9473(03)00007-0](https://doi.org/10.1016/S0167-9473%2803%2900007-0)
.

Bonett DG (2021). *Statistical Methods for Psychologists, Vol 1-5,
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
f <- c(26, 24, 38)
n <- c(60, 60, 60)
v <- c(-.5, -.5, 1)
meta.lc.prop1(.05, f, n, v)
#>           Estimate         SE         LL        UL
#> Contrast 0.2119565 0.07602892 0.06294259 0.3609705

# Should return: 
#           Estimate         SE         LL        UL
# Contrast 0.2119565 0.07602892 0.06294259 0.3609705

```
