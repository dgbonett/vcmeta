# Confidence interval for a log-linear contrast of proportion ratios from 2-group studies

Computes the estimate, standard error, and confidence interval for an
exponentiated log-linear contrast of 2-group proportion ratios from two
or more studies.

For more details, see Section 3.2 of Bonett (2021, Volume 5).

## Usage

``` r
meta.lc.propratio2(alpha, f1, f2, n1, n2, v)
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

- Estimate - estimated log-linear contrast

- SE - standard error of log-linear contrast

- exp(Estimate) - exponentiated log-linear contrast

- exp(LL) - lower limit of the exponentiated confidence interval

- exp(UL) - upper limit of the exponentiated confidence interval

## References

Price RM, Bonett DG (2008). “Confidence intervals for a ratio of two
independent binomial proportions.” *Statistics in Medicine*, **27**(26),
5497–5508. ISSN 02776715,
[doi:10.1002/sim.3376](https://doi.org/10.1002/sim.3376) .

Bonett DG, Price RM (2015). “Varying coefficient meta-analysis methods
for odds ratios and risk ratios.” *Psychological Methods*, **20**(3),
394–406. ISSN 1939-1463,
[doi:10.1037/met0000032](https://doi.org/10.1037/met0000032) .

Bonett DG (2021). *Statistical Methods for Psychologists, Vol 1-5,
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
n1 <- c(50, 150, 150)
f1 <- c(16, 50, 25)
n2 <- c(50, 150, 150)
f2 <- c(7, 15, 20)
v <- c(1, -1, 0)
meta.lc.propratio2(.05, f1, f2, n1, n2, v)
#>            Estimate        SE exp(Estimate)   exp(LL)  exp(UL)
#> Contrast -0.3853396 0.4828218     0.6802196 0.2640405 1.752378

# Should return:
#            Estimate        SE  exp(Estimate)   exp(LL)  exp(UL)
# Contrast -0.3853396 0.4828218      0.6802196 0.2640405 1.752378

```
