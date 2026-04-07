# Confidence interval for a subgroup difference in average effect size

Computes the estimate, standard error, and confidence interval for a
difference in the average effect size (any type of effect size) for two
mutually exclusive subgroups of studies. Each subgroup can have one or
more studies. All of the effects sizes should be compatible.

For more details, see Section 3.3 of Bonett (2021, Volume 5).

## Usage

``` r
meta.sub.gen(alpha, est, se, group)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- est:

  vector of estimated effect sizes

- se:

  vector of effect size standard errors

- group:

  vector of group indicators:

  - 1 for set A

  - 2 for set B

  - 0 to ignore

## Value

Returns a matrix with three rows:

- Row 1 - estimate for Set A

- Row 2 - estimate for Set B

- Row 3 - estimate for difference, Set A - Set B

The columns are:

- Estimate - estimated average effect size or difference

- SE - standard error

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## References

Bonett DG (2021). *Statistical Methods for Psychologists, Vol 1-5,
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
est <- c(.920, .896, .760, .745)
se <- c(.098, .075, .069, .055) 
group <- c(1, 1, 2, 2)
meta.sub.gen(.05, est, se, group)
#>                Estimate         SE          LL        UL
#> Set A:           0.9080 0.06170292 0.787064504 1.0289355
#> Set B:           0.7525 0.04411916 0.666028042 0.8389720
#> Set A - Set B:   0.1555 0.07585348 0.006829917 0.3041701

# Should return:
#                Estimate         SE          LL        UL
# Set A:           0.9080 0.06170292 0.787064504 1.0289355
# Set B:           0.7525 0.04411916 0.666028042 0.8389720
# Set A - Set B:   0.1555 0.07585348 0.006829917 0.3041701

```
