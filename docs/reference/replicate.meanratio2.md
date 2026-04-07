# Compares and combines 2-group mean ratios in original and follow-up studies

This function computes confidence intervals from an original study and a
follow-up study where the effect size is a 2-group mean ratio.
Confidence intervals for the log-difference and average log effect size
are also computed. Equality of variances within or across studies is not
assumed. A Satterthwaite adjustment to the degrees of freedom is used to
improve the accuracy of the confidence intervals. The confidence level
for the difference is 1 – 2\*alpha, which is recommended for equivalence
testing.

For more details, see Chapter 4 of Bonett (2021, Volume 5).

## Usage

``` r
replicate.meanratio2(
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

Bonett DG, Price RM (2020). “Confidence intervals for ratios of means
and medians.” *Journal of Educational and Behavioral Statistics*,
**45**(6), 750–770. ISSN 1076-9986,
[doi:10.3102/1076998620934125](https://doi.org/10.3102/1076998620934125)
.

Bonett DG (2021). *Statistical Methods for Psychologists, Vol 1-5,
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
replicate.meanratio2(.05, 21.9, 16.1, 3.82, 3.21, 40, 40, 25.2, 19.1, 3.98, 3.79, 75, 75)
#>                         Estimate         SE       t       p     df
#> Original:             0.30766736 0.04188600  7.3454 0.00000  76.65
#> Follow-up:            0.27715566 0.02928438  9.4643 0.00000 140.91
#> Original - Follow-up: 0.03051171 0.05110785  0.5970 0.55141 150.35
#> Average:              0.29241151 0.02555393 11.4429 0.00000 150.35
#>                       exp(Estimate)   exp(LL)  exp(UL)
#> Original:                  1.360248 1.2513908 1.478576
#> Follow-up:                 1.319372 1.2451576 1.398009
#> Original - Follow-up:      1.030982 0.9473616 1.121983
#> Average:                   1.339654 1.2736927 1.409032

# Should return:
#                         Estimate         SE       t       p     df
# Original:             0.30766736 0.04188600  7.3454 0.00000  76.65
# Follow-up:            0.27715566 0.02928438  9.4643 0.00000 140.91
# Original - Follow-up: 0.03051171 0.05110785  0.5970 0.55141 150.35
# Average:              0.29241151 0.02555393 11.4429 0.00000 150.35
#                       exp(Estimate)   exp(LL)  exp(UL)
# Original:                  1.360248 1.2513908 1.478576
# Follow-up:                 1.319372 1.2451576 1.398009
# Original - Follow-up:      1.030982 0.9473616 1.121983
# Average:                   1.339654 1.2736927 1.409032

```
