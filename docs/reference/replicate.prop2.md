# Compares and combines 2-group proportion differences in original and follow-up studies

This function computes confidence intervals from an original study and a
follow-up study where the effect size is a 2-group proportion
difference. Confidence intervals for the difference and average effect
size are also computed. The confidence level for the difference is 1 –
2\*alpha, which is recommended for equivalence testing.

For more details, see Chapter 4 of Bonett (2021, Volume 5).

## Usage

``` r
replicate.prop2(alpha, f11, f12, n11, n12, f21, f22, n21, n22)
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

- Row 3 estimates the difference in proportion differences

- Row 4 estimates the average proportion difference

The columns are:

- Estimate - proportion difference estimate (single study, difference,
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
replicate.prop2(.05, 21, 16, 40, 40, 19, 13, 60, 60)
#>                         Estimate         SE     z     p          LL        UL
#> Original:             0.11904762 0.10805233 1.102 0.271 -0.09273105 0.3308263
#> Follow-up:            0.09677419 0.07965047 1.215 0.224 -0.05933787 0.2528863
#> Original - Follow-up: 0.02359056 0.13542107 0.174 0.862 -0.19915727 0.2463384
#> Average:              0.11015594 0.06771053 1.627 0.104 -0.02255427 0.2428661

# Should return:
#                         Estimate         SE     z     p          LL        UL
# Original:             0.11904762 0.10805233 1.102 0.271 -0.09273105 0.3308263
# Follow-up:            0.09677419 0.07965047 1.215 0.224 -0.05933787 0.2528863
# Original - Follow-up: 0.02359056 0.13542107 0.174 0.862 -0.19915727 0.2463384
# Average:              0.11015594 0.06771053 1.627 0.104 -0.02255427 0.2428661

```
