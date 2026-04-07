# Compares and combines 2-group mean differences in original and follow-up studies

This function computes confidence intervals from an original study and a
follow-up study where the effect size is a 2-group mean difference.
Confidence intervals for the difference and average effect size are also
computed. Equality of variances within or across studies is not assumed.
A Satterthwaite adjustment to the degrees of freedom is used to improve
the accuracy of the confidence intervals. The confidence level for the
difference is 1 – 2\*alpha, which is recommended for equivalence
testing.

For more details, see Chapter 4 of Bonett (2021, Volume 5).

## Usage

``` r
replicate.mean2(
  alpha,
  m11,
  m12,
  sd11,
  sd12,
  n11,
  n12,
  m21,
  m22,
  sd21,
  sd22,
  n21,
  n22
)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- m11:

  estimated mean for group 1 in original study

- m12:

  estimated mean for group 2 in original study

- sd11:

  estimated SD for group 1 in original study

- sd12:

  estimated SD for group 2 in original study

- n11:

  sample size for group 1 in original study

- n12:

  sample size for group 2 in original study

- m21:

  estimated mean for group 1 in follow-up study

- m22:

  estimated mean for group 2 in follow-up study

- sd21:

  estimated SD for group 1 in follow-up study

- sd22:

  estimated SD for group 2 in follow-up study

- n21:

  sample size for group 1 in follow-up study

- n22:

  sample size for group 2 in follow-up study

## Value

A 4-row matrix. The rows are:

- Row 1 summarizes the original study

- Row 2 summarizes the follow-up study

- Row 3 estimates the difference in mean differences

- Row 4 estimates the average mean difference

The columns are:

- Estimate - mean difference estimate (single study, difference,
  average)

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
replicate.mean2(.05, 21.9, 16.1, 3.82, 3.21, 40, 40, 
                     25.2, 19.1, 3.98, 3.79, 75, 75)
#>                       Estimate        SE      t     p        LL       UL     df
#> Original:                 5.80 0.7889312  7.352 0.000  4.228624 7.371376  75.75
#> Follow-up:                6.10 0.6346075  9.612 0.000  4.845913 7.354087 147.65
#> Original - Follow-up:    -0.30 1.0124916 -0.296 0.767 -1.974571 1.374571 169.16
#> Average:                  5.95 0.5062458 11.753 0.000  4.950627 6.949373 169.16

# Should return:
#                       Estimate        SE      t     p        LL       UL     df
# Original:                 5.80 0.7889312  7.352 0.000  4.228624 7.371376  75.75
# Follow-up:                6.10 0.6346075  9.612 0.000  4.845913 7.354087 147.65
# Original - Follow-up:    -0.30 1.0124916 -0.296 0.767 -1.974571 1.374571 169.16
# Average:                  5.95 0.5062458 11.753 0.000  4.950627 6.949373 169.16

```
