# Compares and combines Cronbach reliablity in original and follow-up studies

This function can be used to compare and combine a Cronbach reliablity
coefficient from an original study and a follow-up study. The confidence
level for the difference is 1 – 2\*alpha, which is recommended for
equivalence testing.

For more details, see Chapter 4 of Bonett (2021, Volume 5).

## Usage

``` r
replicate.cronbach(alpha, rel1, n1, rel2, n2, r)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- rel1:

  estimated reliability in original study

- n1:

  sample size in original study

- rel2:

  estimated reliability in follow-up study

- n2:

  sample size in follow-up study

- r:

  number of measurements (e.g., items)

## Value

A 4-row matrix. The rows are:

- Row 1 summarizes the original study

- Row 2 summarizes the follow-up study

- Row 3 estimates the difference in correlations

- Row 4 estimates the average correlation

The columns are:

- Estimate - correlation estimate (single study, difference, average)

- SE - standard error

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## References

Bonett DG (2010). “Varying coefficient meta-analytic methods for alpha
reliability.” *Psychological Methods*, **15**(4), 368–385. ISSN
1939-1463, [doi:10.1037/a0020142](https://doi.org/10.1037/a0020142) .

Bonett DG, Price RM (2015). “Varying coefficient meta-analysis methods
for odds ratios and risk ratios.” *Psychological Methods*, **20**(3),
394–406. ISSN 1939-1463,
[doi:10.1037/met0000032](https://doi.org/10.1037/met0000032) .

Bonett DG (2021). “Design and analysis of replication studies.”
*Organizational Research Methods*, **24**(3), 513–529. ISSN 1094-4281,
[doi:10.1177/1094428120911088](https://doi.org/10.1177/1094428120911088)
.

Bonett DG (2021). *Statistical Methods for Psychologists, Vol 1-5,
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
replicate.cronbach(.05, .883, 100, .869, 200, 6)
#>                       Estimate      SE      LL     UL
#> Original:                0.883 0.01831  0.8436 0.9152
#> Follow-up:               0.869 0.01442  0.8387 0.8952
#> Original - Follow-up:    0.014 0.02331 -0.0334 0.0582
#> Average:                 0.876 0.01172  0.8519 0.8977

# Should return:
#                       Estimate      SE      LL     UL
# Original:                0.883 0.01831  0.8436 0.9152
# Follow-up:               0.869 0.01442  0.8387 0.8952
# Original - Follow-up:    0.014 0.02331 -0.0334 0.0582
# Average:                 0.876 0.01172  0.8519 0.8977

```
