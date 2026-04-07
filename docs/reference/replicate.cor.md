# Compares and combines Pearson or partial correlations in original and follow-up studies

This function can be used to compare and combine Pearson or partial
correlations from an original study and a follow-up study. The
confidence level for the difference is 1 – 2\*alpha, which is
recommended for equivalence testing.

For more details, see Chapter 4 of Bonett (2021, Volume 5).

## Usage

``` r
replicate.cor(alpha, cor1, n1, cor2, n2, s)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- cor1:

  estimated correlation in original study

- n1:

  sample size in original study

- cor2:

  estimated correlation in follow-up study

- n2:

  sample size in follow-up study

- s:

  number of control variables in each study (0 for Pearson)

## Value

A 4-row matrix. The rows are:

- Row 1 summarizes the original study

- Row 2 summarizes the follow-up study

- Row 3 estimates the difference in correlations

- Row 4 estimates the average correlation

The columns are:

- Estimate -correlation estimate (single study, difference, average)

- SE - standard error

- z - t-value for rows 1 and 2; z-value for rows 3 and 4

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
replicate.cor(.05, .598, 80, .324, 200, 0)
#>                       Estimate      SE     z     p     LL     UL
#> Original:                0.598 0.07321 6.589 0.000 0.4355 0.7228
#> Follow-up:               0.324 0.06377 4.819 0.000 0.1940 0.4428
#> Original - Follow-up:    0.274 0.09709 2.633 0.008 0.1065 0.4265
#> Average:                 0.461 0.04854 7.635 0.000 0.3725 0.5412

# Should return:
#                       Estimate      SE     z     p     LL     UL
# Original:                0.598 0.07321 6.589 0.000 0.4355 0.7228
# Follow-up:               0.324 0.06377 4.819 0.000 0.1940 0.4428
# Original - Follow-up:    0.274 0.09709 2.633 0.008 0.1065 0.4265
# Average:                 0.461 0.04854 7.635 0.000 0.3725 0.5412

```
