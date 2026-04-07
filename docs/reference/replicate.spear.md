# Compares and combines Spearman correlations in original and follow-up studies

This function can be used to compare and combine Spearman correlations
from an original study and a follow-up study. The confidence level for
the difference is 1 – 2\*alpha, which is recommended for equivalence
testing.

For more details, see Chapter 4 of Bonett (2021, Volume 5).

## Usage

``` r
replicate.spear(alpha, cor1, n1, cor2, n2)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- cor1:

  estimated Spearman correlation in original study

- n1:

  sample size in original study

- cor2:

  estimated Spearman correlation in follow-up study

- n2:

  sample size in follow-up study

## Value

A 4-row matrix. The rows are:

- Row 1 summarizes the original study

- Row 2 summarizes the follow-up study

- Row 3 estimates the difference in correlations

- Row 4 estimates the average correlation

The columns are:

- Estimate - Spearman correlation estimate (single study, difference,
  average)

- SE - standard error

- z - z-value

- p - p-value

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
replicate.spear(.05, .598, 80, .324, 200)
#>                       Estimate      SE     z     p     LL     UL
#> Original:                0.598 0.07948 5.315 0.000 0.4199 0.7318
#> Follow-up:               0.324 0.06542 4.571 0.000 0.1905 0.4457
#> Original - Follow-up:    0.274 0.10294 3.438 0.001 0.0948 0.4342
#> Average:                 0.461 0.05147 9.968 0.000 0.3670 0.5457

# Should return:
#                       Estimate      SE     z     p     LL     UL
# Original:                0.598 0.07948 5.315 0.000 0.4199 0.7318
# Follow-up:               0.324 0.06542 4.571 0.000 0.1905 0.4457
# Original - Follow-up:    0.274 0.10294 3.438 0.001 0.0948 0.4342
# Average:                 0.461 0.05147 9.968 0.000 0.3670 0.5457

```
