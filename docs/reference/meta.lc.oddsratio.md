# Confidence interval for a log-linear contrast of odds ratios

Computes the estimate, standard error, and confidence interval for an
exponentiated log-linear contrast of odds ratios from two or more
studies.

For more details, see Section 3.2 of Bonett (2021, Volume 5).

## Usage

``` r
meta.lc.oddsratio(alpha, f1, f2, n1, n2, v)
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
meta.lc.oddsratio(.05, f1, f2, n1, n2, v)
#>            Estimate        SE exp(Estimate)   exp(LL)  exp(UL)
#> Contrast -0.4596883 0.5895438     0.6314805 0.1988563 2.005305

# Should return:
#            Estimate        SE  exp(Estimate)   exp(LL)  exp(UL)
# Contrast -0.4596883 0.5895438      0.6314805 0.1988563 2.005305

```
