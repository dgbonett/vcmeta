# Compares and combines paired-samples proportion differences in original and follow-up studies

This function computes confidence intervals from an original study and a
follow-up study where the effect size is a paired-samples proportion
difference. Confidence intervals for the difference and average of
effect sizes are also computed. The confidence level for the difference
is 1 – 2\*alpha, which is recommended for equivalence testing.

For more details, see Chapter 4 of Bonett (2021, Volume 5).

## Usage

``` r
replicate.prop.ps(alpha, f1, f2)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- f1:

  vector of frequency counts for 2x2 table in original study

- f2:

  vector of frequency counts for 2x2 table in follow-up study

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

Bonett DG, Price RM (2012). “Adjusted wald confidence interval for a
difference of binomial proportions based on paired data.” *Journal of
Educational and Behavioral Statistics*, **37**(4), 479–488. ISSN
1076-9986,
[doi:10.3102/1076998611411915](https://doi.org/10.3102/1076998611411915)
.

Bonett DG (2021). *Statistical Methods for Psychologists, Vol 1-5,
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
f1 <- c(42, 2, 15, 61)
f2 <- c(69, 5, 31, 145)
replicate.prop.ps(.05, f1, f2)
#>                          Estimate         SE     z     p          LL         UL
#> Original:             0.106557377 0.03440159 3.097 0.002  0.03913151 0.17398325
#> Follow-up:            0.103174603 0.02358274 4.375 0.000  0.05695329 0.14939592
#> Original - Follow-up: 0.003852359 0.04097037 0.094 0.925 -0.06353791 0.07124263
#> Average:              0.105511837 0.02048519 5.151 0.000  0.06536161 0.14566206

# Should return:
#                          Estimate         SE     z     p          LL         UL
# Original:             0.106557377 0.03440159 3.097 0.002  0.03913151 0.17398325
# Follow-up:            0.103174603 0.02358274 4.375 0.000  0.05695329 0.14939592
# Original - Follow-up: 0.003852359 0.04097037 0.094 0.925 -0.06353791 0.07124263
# Average:              0.105511837 0.02048519 5.151 0.000  0.06536161 0.14566206

```
