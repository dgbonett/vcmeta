# Compares and combines odds ratios in original and follow-up studies

This function computes confidence intervals for an odds ratio from an
original study and a follow-up study. Confidence intervals for the ratio
of odds ratios and geometric average odds ratio are also computed. The
confidence level for the ratio of ratios is 1 – 2\*alpha, which is
recommended for equivalence testing.

For more details, see Chapter 4 of Bonett (2021, Volume 5).

## Usage

``` r
replicate.oddsratio(alpha, est1, se1, est2, se2)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- est1:

  estimate of log odds ratio in original study

- se1:

  standard error of log odds ratio in original study

- est2:

  estimate of log odds ratio in follow-up study

- se2:

  standard error of log odds ratio in follow-up study

## Value

A 4-row matrix. The rows are:

- Row 1 summarizes the original study

- Row 2 summarizes the follow-up study

- Row 3 estimates the ratio of odds ratios

- Row 4 estimates the geometric average odds ratio

The columns are:

- Estimate - log odds ratio estimate (single study, ratio, average)

- SE - standard error of log odds estimate

- z - z-value

- p - p-value

- exp(Estimate) - exponentiated estimate

- exp(LL) - exponentiated lower limit of the confidence interval

- exp(UL) - exponentiated upper limit of the confidence interval

## References

Bonett DG (2021). “Design and analysis of replication studies.”
*Organizational Research Methods*, **24**(3), 513–529. ISSN 1094-4281,
[doi:10.1177/1094428120911088](https://doi.org/10.1177/1094428120911088)
.

Bonett DG (2021). *Statistical Methods for Psychologists, Vol 1-5,
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
replicate.oddsratio(.05, 1.39, .302, 1.48, .206)
#>                        Estimate        SE      z     p exp(Estimate)   exp(LL)
#> Original:            1.39000000 0.3020000  4.603 0.000     4.0148501 2.2212961
#> Follow-up:           1.48000000 0.2060000  7.184 0.000     4.3929457 2.9336501
#> Original/Follow-up: -0.06273834 0.3655681 -0.172 0.864     0.9391892 0.5147653
#> Average:             0.36067292 0.1827840  1.973 0.048     1.4342943 1.0024257
#>                      exp(UL)
#> Original:           7.256583
#> Follow-up:          6.578144
#> Original/Follow-up: 1.713551
#> Average:            2.052222

# Should return:
#                        Estimate        SE      z     p
# Original:            1.39000000 0.3020000  4.603 0.000
# Follow-up:           1.48000000 0.2060000  7.184 0.000
# Original/Follow-up: -0.06273834 0.3655681 -0.172 0.864
# Average:             0.36067292 0.1827840  1.973 0.048
#                     exp(Estimate)   exp(LL)  exp(UL)
# Original:               4.0148501 2.2212961 7.256583
# Follow-up:              4.3929457 2.9336501 6.578144
# Original/Follow-up:     0.9391892 0.5147653 1.713551
# Average:                1.4342943 1.0024257 2.052222

```
