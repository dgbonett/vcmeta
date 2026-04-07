# Compares and combines single proportions in original and follow-up studies

This function computes confidence intervals for a single proportion from
an original study and a follow-up study. Confidence intervals for the
difference between the two proportions and average of the two
proportions are also computed. The confidence level for the difference
is 1 – 2\*alpha, which is recommended for equivalence testing.

For more details, see Chapter 4 of Bonett (2021, Volume 5).

## Usage

``` r
replicate.prop1(alpha, f1, n1, f2, n2)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- f1:

  frequency count in original study

- n1:

  sample size in original study

- f2:

  frequency count in follow-up study

- n2:

  sample size for in follow-up study

## Value

A 4-row matrix. The rows are:

- Row 1 summarizes the original study

- Row 2 summarizes the follow-up study

- Row 3 estimates the difference in proportions

- Row 4 estimates the average proportion

The columns are:

- Estimate - proportion estimate (single study, difference, average)

- SE - standard error

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## References

Bonett DG (2021). “Design and analysis of replication studies.”
*Organizational Research Methods*, **24**(3), 513–529. ISSN 1094-4281,
[doi:10.1177/1094428120911088](https://doi.org/10.1177/1094428120911088)
.

Bonett DG (2021). *Statistical Methods for Psychologists, Vol 1-5,
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
replicate.prop1(.05, 21, 300, 35, 400)
#>                          Estimate         SE          LL         UL
#> Original:              0.07565789 0.01516725  0.04593064 0.10538515
#> Follow-up:             0.09158416 0.01435033  0.06345803 0.11971029
#> Original - Follow-up: -0.01670456 0.02065098 -0.05067239 0.01726328
#> Average:               0.08119996 0.01032549  0.06096237 0.10143755

# Should return:
#                          Estimate         SE          LL         UL
# Original:              0.07565789 0.01516725  0.04593064 0.10538515
# Follow-up:             0.09158416 0.01435033  0.06345803 0.11971029
# Original - Follow-up: -0.01670456 0.02065098 -0.05067239 0.01726328
# Average:               0.08119996 0.01032549  0.06096237 0.10143755

```
