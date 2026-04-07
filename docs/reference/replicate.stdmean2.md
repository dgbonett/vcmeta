# Compares and combines 2-group standardized mean differences in original and follow-up studies

This function computes confidence intervals from an original study and a
follow-up study where the effect size is a 2-group standardized mean
difference. Confidence intervals for the difference and average effect
size are also computed. Equality of variances within or across studies
is not assumed. The confidence level for the difference is 1 – 2\*alpha,
which is recommended for equivalence testing. Square root unweighted
variances, square root weighted variances, and single-group standard
deviation are options for the standardizer.

For more details, see Chapter 4 of Bonett (2021, Volume 5).

## Usage

``` r
replicate.stdmean2(
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
  n22,
  stdzr
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

- stdzr:

  - set to 0 for square root unweighted average variance standardizer

  - set to 1 for group 1 SD standardizer

  - set to 2 for group 2 SD standardizer

  - set to 3 for square root weighted average variance standardizer

## Value

A 4-row matrix. The rows are:

- Row 1 summarizes the original study

- Row 2 summarizes the follow-up study

- Row 3 estimates the difference in standardized mean differences

- Row 4 estimates the average standardized mean difference

The columns are:

- Estimate - standardized mean difference estimate (single study,
  difference, average)

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
replicate.stdmean2(.05, 21.9, 16.1, 3.82, 3.21, 40, 40, 
                        25.2, 19.1, 3.98, 3.79, 75, 75, 0)
#>                       Estimate     SE      LL     UL
#> Original:               1.6280 0.2629  1.1286 2.1592
#> Follow-up:              1.5617 0.1881  1.2011 1.9383
#> Original - Follow-up:   0.0742 0.3232 -0.4575 0.6059
#> Average:                1.5949 0.1616  1.2781 1.9116

# Should return: 
#                        Estimate     SE      LL     UL
#  Original:               1.6280 0.2595  1.1353 2.1524
#  Follow-up:              1.5617 0.1871  1.2030 1.9363
#  Original - Follow-up:   0.0742 0.3199 -0.4519 0.6004
#  Average:                1.5949 0.1599  1.2814 1.9083

```
