# Compares and combines point-biserial correlations in original and follow-up studies

This function computes confidence intervals from an original study and a
follow-up study where the effect size is a point-biserial correlation.
Confidence intervals for the difference and average effect size are also
computed. The confidence level for the difference is 1 – 2\*alpha, which
is recommended for equivalence testing. The point-biserial correlation
in each study is computed from a standardized mean difference. Two types
of standardized mean differences can be requested. One type uses the
square root of unweighted variances as a standardizer and is recommended
for 2-group experimental designs. The other type uses the square root of
weighted variances as a standardizer and is recommended for 2-group
non-experimental designs with simple random sampling. Equality of
variances across or within studies is not assumed.

For more details, see Chapter 4 of Bonett (2021, Volume 5).

## Usage

``` r
replicate.pbcor(
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
  type
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

- type:

  - set to 1 for square root weighted average variance standardizer

  - set to 2 for square root unweighted average variance standardizer

## Value

A 4-row matrix. The rows are:

- Row 1 summarizes the original study

- Row 2 summarizes the follow-up study

- Row 3 estimates the difference in point-biserial correlations

- Row 4 estimates the average point-biserial correlation

The columns are:

- Estimate - point-biserial correlation (single study, difference,
  average)

- SE - standard error

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## References

Bonett DG (2020). “Point-biserial correlation: Interval estimation,
hypothesis testing, meta-analysis, and sample size determination.”
*British Journal of Mathematical and Statistical Psychology*,
**73**(S1), 113–144. ISSN 0007-1102,
[doi:10.1111/bmsp.12189](https://doi.org/10.1111/bmsp.12189) .

Bonett DG (2021). *Statistical Methods for Psychologists, Vol 1-5,
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
replicate.pbcor(.05, 21.9, 16.1, 3.82, 3.21, 40, 40, 25.2, 19.1, 3.98, 3.79, 75, 75, 2)
#>                       Estimate      SE      LL     UL
#> Original:               0.6350 0.06061  0.4915 0.7336
#> Follow-up:              0.6174 0.04578  0.5148 0.6959
#> Original - Follow-up:   0.0176 0.07595 -0.1460 0.1599
#> Average:                0.6262 0.03798  0.5460 0.6950

# Should return: 
#                       Estimate      SE      LL     UL
# Original:               0.6350 0.06061  0.4915 0.7336
# Follow-up:              0.6174 0.04578  0.5148 0.6959
# Original - Follow-up:   0.0176 0.07595 -0.1460 0.1599
# Average:                0.6262 0.03798  0.5460 0.6950

replicate.pbcor(.05, 12.2, 10.4, 1.74, 1.59, 68, 94, 13.0, 10.9, 1.48, 1.29, 124, 189, 1)
#>                       Estimate      SE      LL      UL
#> Original:               0.4753 0.05847  0.3487  0.5781
#> Follow-up:              0.6016 0.03365  0.5292  0.6617
#> Original - Follow-up:  -0.1262 0.06746 -0.2664 -0.0005
#> Average:                0.5384 0.03373  0.4691  0.6012

# Should return: 
#                       Estimate      SE      LL      UL
# Original:               0.4753 0.05847  0.3487  0.5781
# Follow-up:              0.6016 0.03365  0.5292  0.6617
# Original - Follow-up:  -0.1262 0.06746 -0.2664 -0.0005
# Average:                0.5384 0.03373  0.4691  0.6012

```
