# Compares and combines paired-samples mean ratios in original and follow-up studies

This function computes confidence intervals from an original study and a
follow-up study where the effect size is a paired-samples mean ratio.
Confidence intervals for the log-difference and average log effect size
are also computed. Equality of variances within or across studies is not
assumed. A Satterthwaite adjustment to the degrees of freedom is used to
improve the accuracy of the confidence intervals. The confidence level
for the difference is 1 – 2\*alpha, which is recommended for equivalence
testing.

For more details, see Chapter 4 of Bonett (2021, Volume 5).

## Usage

``` r
replicate.meanratio.ps(
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

  estimated mean for group 1 in original study

- m12:

  estimated mean for group 2 in original study

- sd11:

  estimated SD for group 1 in original study

- sd12:

  estimated SD for group 2 in original study

- cor1:

  estimated correlation of paired measurements in original study

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
replicate.meanratio.ps(.05, 86.22, 70.93, 14.89, 12.32, .765, 20, 
                       84.81, 77.24, 15.68, 16.95, .702, 75)
#>                        Estimate         SE      t       p    df exp(Estimate)
#> Original:             0.1952087 0.02655112 7.3522 0.00000 19.00      1.215565
#> Follow-up:            0.0934960 0.01839400 5.0830 0.00000 74.00      1.098006
#> Original - Follow-up: 0.1017127 0.03230018 3.1490 0.00313 39.29      1.107065
#> Average:              0.1443523 0.01615009 8.9382 0.00000 39.29      1.155291
#>                        exp(LL)  exp(UL)
#> Original:             1.149856 1.285028
#> Follow-up:            1.058492 1.138996
#> Original - Follow-up: 1.048437 1.168972
#> Average:              1.118170 1.193645

# Should return:
#                         Estimate         SE       t       p    df
# Original:             0.30766736 0.02420564 12.7106 0.00000 39.00
# Follow-up:            0.27715566 0.01532870 18.0808 0.00000 74.00
# Original - Follow-up: 0.03051171 0.02865104  1.0649 0.29055 70.57
# Average:              0.29241151 0.01432552 20.4119 0.00000 70.57
#                       exp(Estimate)   exp(LL)  exp(UL)
# Original:                  1.360248 1.2952540 1.428504
# Follow-up:                 1.319372 1.2796832 1.360291
# Original - Follow-up:      1.030982 0.9829058 1.081410
# Average:                   1.339654 1.3019254 1.378476

```
