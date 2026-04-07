# Compares and combines any type of correlation in original and follow-up studies

This function can be used to compare and combine any type of correlation
from an original study and a follow-up study. The confidence level for
the difference is 1 – 2\*alpha, which is recommended for equivalence
testing.

For more details, see Chapter 4 of Bonett (2021, Volume 5).

## Usage

``` r
replicate.cor.gen(alpha, cor1, se1, cor2, se2)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- cor1:

  estimated correlation in original study

- se1:

  standard error of correlation in original study

- cor2:

  estimated correlation in follow-up study

- se2:

  standard error of correlation in follow-up study

## Value

A 4-row matrix. The rows are:

- Row 1 summarizes the original study

- Row 2 summarizes the follow-up study

- Row 3 estimates the difference in correlations

- Row 4 estimates the average correlation

The columns are:

- Estimate - correlation estimate (single study, difference, average)

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
replicate.cor.gen(.05, .454, .170, .318, .098)
#>                       Estimate      SE     z     p      LL     UL
#> Original:                0.454 0.17000 2.287 0.022  0.0699 0.7209
#> Follow-up:               0.318 0.09800 3.022 0.003  0.1152 0.4953
#> Original - Follow-up:    0.136 0.19622 0.667 0.505 -0.2154 0.4237
#> Average:                 0.386 0.09811 3.409 0.001  0.1961 0.5480

# Should return:
#                       Estimate      SE     z     p      LL     UL
# Original:                0.454 0.17000 2.287 0.022  0.0699 0.7209
# Follow-up:               0.318 0.09800 3.022 0.003  0.1152 0.4953
# Original - Follow-up:    0.136 0.19622 0.667 0.505 -0.2154 0.4237
# Average:                 0.386 0.09811 3.409 0.001  0.1961 0.5480

```
