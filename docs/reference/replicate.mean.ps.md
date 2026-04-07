# Compares and combines paired-samples mean differences in original and follow-up studies

This function computes confidence intervals from an original study and a
follow-up study where the effect size is a paired-samples mean
difference. Confidence intervals for the difference and average effect
size are also computed. Equality of variances within or across studies
is not assumed. A Satterthwaite adjustment to the degrees of freedom is
used to improve the accuracy of the confidence intervals for the
difference and average. The confidence level for the difference is 1 –
2\*alpha, which is recommended for equivalence testing.

For more details, see Chapter 4 of Bonett (2021, Volume 5).

## Usage

``` r
replicate.mean.ps(
  alpha,
  m11,
  m12,
  sd11,
  sd12,
  cor1,
  n1,
  m21,
  m22,
  sd21,
  sd22,
  cor2,
  n2
)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- m11:

  estimated mean for measurement 1 in original study

- m12:

  estimated mean for measurement 2 in original study

- sd11:

  estimated SD for measurement 1 in original study

- sd12:

  estimated SD for measurement 2 in original study

- cor1:

  estimated correlation of paired measurements in orginal study

- n1:

  sample size in original study

- m21:

  estimated mean for measurement 1 in follow-up study

- m22:

  estimated mean for measurement 2 in follow-up study

- sd21:

  estimated SD for measurement 1 in follow-up study

- sd22:

  estimated SD for measurement 2 in follow-up study

- cor2:

  estimated correlation of paired measurements in follow-up study

- n2:

  sample size in follow-up study

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

- df - degrees of freedom

- t - t-value

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
replicate.mean.ps(.05, 86.22, 70.93, 14.89, 12.32, .765, 20, 
                       84.81, 77.24, 15.68, 16.95, .702, 75)
#>                       Estimate       SE     t     p        LL       UL   df
#> Original:                15.29 2.154344 7.097 0.000 10.780906 19.79909 19.0
#> Follow-up:                7.57 1.460664 5.183 0.000  4.659564 10.48044 74.0
#> Original - Follow-up:     7.72 2.602832 2.966 0.005  3.332885 12.10712 38.4
#> Average:                 11.43 1.301416 8.783 0.000  8.796322 14.06368 38.4

#  Should return:
#                        Estimate       SE     t     p        LL       UL   df
#  Original:                15.29 2.154344 7.097 0.000 10.780906 19.79909 19.0
#  Follow-up:                7.57 1.460664 5.183 0.000  4.659564 10.48044 74.0
#  Original - Follow-up:     7.72 2.602832 2.966 0.005  3.332885 12.10712 38.4
#  Average:                 11.43 1.301416 8.783 0.000  8.796322 14.06368 38.4

```
