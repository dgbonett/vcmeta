# Compares and combines paired-samples standardized mean differences in original and follow-up studies

This function computes confidence intervals from an original study and a
follow-up study where the effect size is a paired-samples standardized
mean difference. Confidence intervals for the difference and average
effect size are also computed. Equality of variances within or across
studies is not assumed. The confidence level for the difference is 1 –
2\*alpha, which is recommended for equivalence testing. Square root
unweighted variances and single-condition standard deviation are options
for the standardizer.

For more details, see Chapter 4 of Bonett (2021, Volume 5).

## Usage

``` r
replicate.stdmean.ps(
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
  n2,
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

- cor1:

  estimated correlation of paired observations in orginal study

- n1:

  sample size in original study

- m21:

  estimated mean for group 1 in follow-up study

- m22:

  estimated mean for group 2 in follow-up study

- sd21:

  estimated SD for group 1 in follow-up study

- sd22:

  estimated SD for group 2 in follow-up study

- cor2:

  estimated correlation of paired observations in follow-up study

- n2:

  sample size in follow-up study

- stdzr:

  - set to 0 for square root unweighted average variance standardizer

  - set to 1 for measurement 1 SD standardizer

  - set to 2 for measurement 2 SD standardizer

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
replicate.stdmean.ps(alpha = .05, 86.22, 70.93, 14.89, 12.32, .765, 20, 
                                  84.81, 77.24, 15.68, 16.95, .702, 75, 0)
#>                       Estimate     SE     LL     UL
#> Orginal:                1.0890 0.2292 0.6697 1.5680
#> Follow-up:              0.4605 0.0959 0.2757 0.6516
#> Original - Follow-up:   0.6552 0.2484 0.2466 1.0638
#> Average:                0.7748 0.1242 0.5313 1.0182

# Should return:
#                       Estimate     SE     LL     UL
# Orginal:                1.0890 0.2292 0.6697 1.5680
# Follow-up:              0.4605 0.0959 0.2757 0.6516
# Original - Follow-up:   0.6552 0.2484 0.2466 1.0638
# Average:                0.7748 0.1242 0.5313 1.0182

```
