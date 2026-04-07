# Compares and combines slope coefficients in original and follow-up studies

This function computes confidence intervals for an OLS slope in a GLM
from the original and follow-up studies, the difference in slopes, and
the average of the slopes. Equality of error variances across studies is
not assumed. The confidence interval for the difference uses a 1 -
2\*alpha confidence level, which is recommended for equivalence testing.
Use the
[replicate.gen](https://dgbonett.github.io/vcmeta/reference/replicate.gen.md)
function for slopes in other types of models (e.g., binary logistic,
ordinal logistic, SEM). A Satterthwaite adjustment to the degrees of
freedom is used to improve the accuracy of the confidence intervals for
the average and the difference.

For more details, see Chapter 4 of Bonett (2021, Volume 5).

## Usage

``` r
replicate.slope(alpha, b1, se1, n1, b2, se2, n2, s)
```

## Arguments

- alpha:

  alpha level for 1-alpha or 1 - 2alpha confidence

- b1:

  sample slope in original study

- se1:

  standard error of slope in original study

- n1:

  sample size in original study

- b2:

  sample slope in follow-up study

- se2:

  standard error of slope in follow-up study

- n2:

  sample size in follow-up study

- s:

  number of predictor variables in model

## Value

A 4-row matrix. The rows are:

- Row 1 summarizes the original study

- Row 2 summarizes the follow-up study

- Row 3 estimates the difference in slopes

- Row 4 estimates the average slope

The columns are:

- Estimate - slope estimate (single study, difference, average)

- SE - standard error

- t - t-value

- p - p-value

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

- df - degrees of freedom

## References

Bonett DG (2021). “Design and analysis of replication studies.”
*Organizational Research Methods*, **24**(3), 513–529. ISSN 1094-4281,
[doi:10.1177/1094428120911088](https://doi.org/10.1177/1094428120911088)
.

Bonett DG (2021). *Statistical Methods for Psychologists, Vol 1-5,
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
replicate.slope(.05, 23.4, 5.16, 50, 18.5, 4.48, 90, 4)
#>                       Estimate       SE     t     p        LL       UL    df
#> Original:                23.40 5.160000 4.535 0.000 13.007227 33.79277  45.0
#> Follow-up:               18.50 4.480000 4.129 0.000  9.592560 27.40744  85.0
#> Original - Follow-up:     4.90 6.833447 0.717 0.475 -6.438743 16.23874 106.4
#> Average:                 20.95 3.416724 6.132 0.000 14.176310 27.72369 106.4

# Should return: 
#                       Estimate       SE     t     p        LL       UL    df
# Original:                23.40 5.160000 4.535 0.000 13.007227 33.79277  45.0
# Follow-up:               18.50 4.480000 4.129 0.000  9.592560 27.40744  85.0
# Original - Follow-up:     4.90 6.833447 0.717 0.475 -6.438743 16.23874 106.4
# Average:                 20.95 3.416724 6.132 0.000 14.176310 27.72369 106.4

```
