# Compares and combines effect sizes in original and follow-up studies

This function can be used to compare and combine any effect size using
the effect size estimate and its standard error from the original study
and the follow-up study. The confidence level for the difference is 1 –
2\*alpha, which is recommended for equivalence testing.

For more details, see Chapter 4 of Bonett (2021, Volume 5).

## Usage

``` r
replicate.gen(alpha, est1, se1, est2, se2)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- est1:

  estimated effect size in original study

- se1:

  effect size standard error in original study

- est2:

  estimated effect size in follow-up study

- se2:

  effect size standard error in follow-up study

## Value

A 4-row matrix. The rows are:

- Row 1 summarizes the original study

- Row 2 summarizes the follow-up study

- Row 3 estimates the difference in effect sizes

- Row 4 estimates the average effect size

Columns are:

- Estimate - effect size estimate (single study, difference, average)

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

Bonett DG (2021). *Statistical Methods for Psychologists, Vol 1-5,
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
replicate.gen(.05, .782, .210, .650, .154)
#>                       Estimate        SE     z     p         LL        UL
#> Original:                0.782 0.2100000 3.724 0.000  0.3704076 1.1935924
#> Follow-up:               0.650 0.1540000 4.221 0.000  0.3481655 0.9518345
#> Original - Follow-up:    0.132 0.2604151 0.507 0.612 -0.2963446 0.5603446
#> Average:                 0.716 0.1302075 5.499 0.000  0.4607979 0.9712021

# Should return: 
#                       Estimate        SE     z     p         LL        UL
# Original:                0.782 0.2100000 3.724 0.000  0.3704076 1.1935924
# Follow-up:               0.650 0.1540000 4.221 0.000  0.3481655 0.9518345
# Original - Follow-up:    0.132 0.2604151 0.507 0.612 -0.2963446 0.5603446
# Average:                 0.716 0.1302075 5.499 0.000  0.4607979 0.9712021

```
