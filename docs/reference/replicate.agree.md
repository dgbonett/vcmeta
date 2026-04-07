# Compares and combines G-index of agreement in original and follow-up studies

This function computes adjusted Wald confidence intervals from an
original study and a follow-up study where the effect size is a G-index
of agreement. Adjusted Wald confidence intervals for the difference and
average effect size are also computed. The confidence level for the
difference is 1 – 2\*alpha, which is recommended for equivalence
testing. As a measurement of agreement, the G-index is usually preferred
to Cohen's kappa.

For more details, see Chapter 4 of Bonett (2021, Volume 5).

## Usage

``` r
replicate.agree(alpha, f1, n1, f2, n2, k)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- f1:

  number of objects rated in agreement in original study

- n1:

  sample size (number of objects) in original study

- f2:

  number of objects rated in agreement in follow-up study

- n2:

  sample size (number of objects) in follow-up study

- k:

  number of rating categories

## Value

A 4-row matrix. The rows are:

- Row 1 summarizes the original study

- Row 2 summarizes the follow-up study

- Row 3 estimates the difference in G-indicies

- Row 4 estimates the average G-index

The columns are:

- Estimate - MLE of G-index (single study, difference, average)

- SE - standard error of adjusted estimate

- LL - lower limit of the adjusted confidence interval

- UL - upper limit of the adjusted confidence interval

## References

Bonett DG (2022). “Statistical inference for G-indices of agreement.”
*Journal of Educational and Behavioral Statistics*, **47**(4), 438–458.
[doi:10.3102/10769986221088561](https://doi.org/10.3102/10769986221088561)
.

Bonett DG (2021). “Design and analysis of replication studies.”
*Organizational Research Methods*, **24**(3), 513–529. ISSN 1094-4281,
[doi:10.1177/1094428120911088](https://doi.org/10.1177/1094428120911088)
.

Bonett DG (2021). *Statistical Methods for Psychologists, Vol 1-5,
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
replicate.agree(.05, 85, 100, 160, 200, 2)
#>                       Estimate      SE      LL     UL
#> Original:                 0.70 0.07252  0.5309 0.8152
#> Follow-up:                0.60 0.05662  0.4773 0.6992
#> Original - Follow-up:     0.10 0.09160 -0.0584 0.2429
#> Average:                  0.65 0.04580  0.5504 0.7299

# Should return:
#                       Estimate      SE      LL     UL
# Original:                 0.70 0.07252  0.5309 0.8152
# Follow-up:                0.60 0.05662  0.4773 0.6992
# Original - Follow-up:     0.10 0.09160 -0.0584 0.2429
# Average:                  0.65 0.04580  0.5504 0.7299

```
