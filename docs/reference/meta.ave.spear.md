# Confidence interval for an average Spearman correlation

Computes the estimate, standard error, and confidence interval for an
average Spearman correlation from two or more studies. The Spearman
correlation is preferred to the Pearson correlation if the relation
between the two quantitative variables is monotonic rather than linear
or if the bivariate normality assumption is not plausible.

For more details, see Chapter 2 of Bonett (2021, Volume 5).

## Usage

``` r
meta.ave.spear(alpha, n, cor, bystudy = TRUE)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- n:

  vector of sample sizes

- cor:

  vector of estimated Spearman correlations

- bystudy:

  logical to also return each study estimate (TRUE) or not

## Value

Returns a matrix. The first row is the average estimate across all
studies. If bystudy is TRUE, there is 1 additional row for each study.
The matrix has the following columns:

- Estimate - estimated effect size

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
n <- c(150, 200, 300, 200, 350)
cor <- c(.14, .29, .16, .21, .23)
meta.ave.spear(.05, n, cor, bystudy = TRUE)
#>         Estimate      SE      LL     UL
#> Average    0.206 0.02944  0.1476 0.2629
#> Study 1    0.140 0.08071 -0.0215 0.2944
#> Study 2    0.290 0.06628  0.1548 0.4146
#> Study 3    0.160 0.05671  0.0469 0.2691
#> Study 4    0.210 0.06850  0.0719 0.3402
#> Study 5    0.230 0.05136  0.1269 0.3282

# Should return:
#         Estimate      SE      LL     UL
# Average    0.206 0.02944  0.1476 0.2629
# Study 1    0.140 0.08071 -0.0215 0.2944
# Study 2    0.290 0.06628  0.1548 0.4146
# Study 3    0.160 0.05671  0.0469 0.2691
# Study 4    0.210 0.06850  0.0719 0.3402
# Study 5    0.230 0.05136  0.1269 0.3282

```
