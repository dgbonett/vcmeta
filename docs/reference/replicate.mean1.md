# Compares and combines single mean in original and follow-up studies

This function computes confidence intervals for a single mean from an
original study and a follow-up study. Confidence intervals for the
difference between the two means and average of the two means are also
computed. Equality of variances across studies is not assumed. A
Satterthwaite adjustment to the degrees of freedom is used to improve
the accuracy of the confidence intervals for the difference and average.
The confidence level for the difference is 1 – 2\*alpha, which is
recommended for equivalence testing.

For more details, see Chapter 4 of Bonett (2021, Volume 5).

## Usage

``` r
replicate.mean1(alpha, m1, sd1, n1, m2, sd2, n2)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- m1:

  estimated mean in original study

- sd1:

  estimated SD in original study

- n1:

  sample size in original study

- m2:

  estimated mean in follow-up study

- sd2:

  estimated SD in follow-up study

- n2:

  sample size for in follow-up study

## Value

A 4-row matrix. The rows are:

- Row 1 summarizes the original study

- Row 2 summarizes the follow-up study

- Row 3 estimates the difference in means

- Row 4 estimates the average mean

The columns are:

- Estimate - mean estimate (single study, difference, average)

- SE - standard error

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
replicate.mean1(.05, 21.9, 3.82, 40, 25.2, 3.98, 75)
#>                       Estimate        SE        LL        UL     df
#> Original:                21.90 0.6039950 20.678305 23.121695 39.000
#> Follow-up:               25.20 0.4595708 24.284285 26.115715 74.000
#> Original - Follow-up:    -3.30 0.7589567 -4.562527 -2.037473 82.630
#> Average:                 23.55 0.3794784 22.795183 24.304817 82.633

# Should return:
#                       Estimate        SE        LL        UL    df
# Original:                21.90 0.6039950 20.678305 23.121695 39.00
# Follow-up:               25.20 0.4595708 24.284285 26.115715 74.00
# Original - Follow-up:    -3.30 0.7589567 -4.562527 -2.037473 82.63
# Average:                 23.55 0.3794784 22.795183 24.304817 82.63

```
