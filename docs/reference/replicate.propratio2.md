# Compares and combines 2-group proportion ratios in original and follow-up studies

This function computes confidence intervals from an original study and a
follow-up study where the effect size is a 2-group log proportion ratio.
Confidence intervals for exponentiated effect sizes are also computed.
The confidence level for the ratio of ratios is 1 – 2\*alpha, which is
recommended for equivalence testing.

For more details, see Chapter 4 of Bonett (2021, Volume 5).

## Usage

``` r
replicate.propratio2(alpha, f11, f12, n11, n12, f21, f22, n21, n22)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- f11:

  frequency count for group 1 in original study

- f12:

  frequency count for group 2 in original study

- n11:

  sample size for group 1 in original study

- n12:

  sample size for group 2 in original study

- f21:

  frequency count for group 1 in follow-up study

- f22:

  frequency count for group 2 in follow-up study

- n21:

  sample size for group 1 in follow-up study

- n22:

  sample size for group 2 in follow-up study

## Value

A 4-row matrix. The rows are:

- Row 1 summarizes the original study

- Row 2 summarizes the follow-up study

- Row 3 estimates the ratio of proportion ratios

- Row 4 estimates the geometric average proportion ratio

The columns are:

- Estimate - log proportion ratio estimate (single study, ratio,
  average)

- SE - standard error of log proportion estimate

- z - z-value

- p - p-value

- exp(Estimate) - exponentiated estimate

- exp(LL) - exponentiated lower limit of the confidence interval

- exp(UL) - exponentiated upper limit of the confidence interval

## References

Bonett DG (2021). *Statistical Methods for Psychologists, Vol 1-5,
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
replicate.propratio2(.05, 21, 16, 40, 40, 19, 13, 60, 60)
#>                       Estimate        SE      z     p exp(Estimate)   exp(LL)
#> Original:            0.2682640 0.2463597  1.089 0.276     1.3076923 0.8068705
#> Follow-up:           0.3735135 0.3082711  1.212 0.226     1.4528302 0.7939881
#> Original/Follow-up: -0.1052495 0.3946190 -0.267 0.789     0.9000999 0.4703209
#> Average:             0.3208887 0.1973095  1.626 0.104     1.3783522 0.9362893
#>                      exp(UL)
#> Original:           2.119373
#> Follow-up:          2.658372
#> Original/Follow-up: 1.722611
#> Average:            2.029132

# Should return:
#                       Estimate        SE      z     p
# Original:            0.2682640 0.2463597  1.089 0.276
# Follow-up:           0.3735135 0.3082711  1.212 0.226
# Original/Follow-up: -0.1052495 0.3946190 -0.267 0.789
# Average:             0.3208887 0.1973095  1.626 0.104
#                     exp(Estimate)   exp(LL)  exp(UL)
# Original:               1.3076923 0.8068705 2.119373
# Follow-up:              1.4528302 0.7939881 2.658372
# Original/Follow-up:     0.9000999 0.4703209 1.722611
# Average:                1.3783522 0.9362893 2.029132

```
