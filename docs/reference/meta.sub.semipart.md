# Confidence interval for a subgroup difference in average semipartial correlations

Computes the estimate, standard error, and confidence interval for a
difference in average semipartial correlations for two subgroups of
mutually exclusive studies. Each subgroup can have one or more studies.

For more details, see Section 3.3 of Bonett (2021, Volume 5).

## Usage

``` r
meta.sub.semipart(alpha, n, cor, r2, group)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- n:

  vector of sample sizes

- cor:

  vector of estimated semi-partial correlations

- r2:

  vector of squared multiple correlations for a model that includes the
  IV and all control variables

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

- Estimate - estimated average correlation or difference

- SE - standard error

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## References

Bonett DG (2021). *Statistical Methods for Psychologists, Vol 1-5,
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
n <- c(55, 190, 65, 35)
cor <- c(.40, .65, .60, .45)
r2 <- c(.25, .41, .43, .39)
group <- c(1, 1, 2, 0)  
meta.sub.semipart(.05, n, cor, r2, group)
#>                Estimate      SE      LL     UL
#> Set A:            0.525 0.05955  0.3987 0.6318
#> Set B:            0.600 0.07931  0.4221 0.7334
#> Set A - Set B:   -0.075 0.09918 -0.2587 0.1325

# Should return:
#                Estimate      SE      LL     UL
# Set A:            0.525 0.05955  0.3987 0.6318
# Set B:            0.600 0.07931  0.4221 0.7334
# Set A - Set B:   -0.075 0.09918 -0.2587 0.1325

```
