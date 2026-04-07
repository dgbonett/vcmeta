# Confidence interval for a subgroup difference in average Spearman correlations

Computes the estimate, standard error, and confidence interval for a
difference in average Spearman correlations for two mutually exclusive
subgroups of studies. Each subgroup can have one or more studies.

For more details, see Section 3.3 of Bonett (2021, Volume 5).

## Usage

``` r
meta.sub.spear(alpha, n, cor, group)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- n:

  vector of sample sizes

- cor:

  vector of estimated Spearman correlations

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

Bonett DG (2008). “Meta-analytic interval estimation for bivariate
correlations.” *Psychological Methods*, **13**(3), 173–181. ISSN
1939-1463, [doi:10.1037/a0012868](https://doi.org/10.1037/a0012868) .

Bonett DG (2021). *Statistical Methods for Psychologists, Vol 1-5,
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
n <- c(55, 190, 65, 35)
cor <- c(.40, .65, .60, .45)
group <- c(1, 1, 2, 0)
meta.sub.spear(.05, n, cor, group)
#>                Estimate      SE      LL     UL
#> Set A:            0.525 0.06484  0.3866 0.6403
#> Set B:            0.600 0.08829  0.3992 0.7459
#> Set A - Set B:   -0.075 0.10954 -0.2761 0.1565

# Should return:
#                Estimate      SE      LL     UL
# Set A:            0.525 0.06484  0.3866 0.6403
# Set B:            0.600 0.08829  0.3992 0.7459
# Set A - Set B:   -0.075 0.10954 -0.2761 0.1565

```
