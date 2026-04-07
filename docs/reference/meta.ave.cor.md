# Confidence interval for an average Pearson or partial correlation

Computes the estimate, standard error, and confidence interval for an
average Pearson or partial correlation from two or more studies. The
sample correlations must be all Pearson correlations or all partial
correlations. Use the
[meta.ave.cor.gen](https://dgbonett.github.io/vcmeta/reference/meta.ave.cor.gen.md)
function to meta-analyze any combination of Pearson, partial, or
Spearman correlations.

For more details, see Chapter 2 of Bonett (2021, Volume 5).

## Usage

``` r
meta.ave.cor(alpha, n, cor, s, bystudy = TRUE)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- n:

  vector of sample sizes

- cor:

  vector of estimated correlations

- s:

  number of control variables (set to 0 for Pearson)

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
n <- c(55, 190, 65, 35)
cor <- c(.40, .65, .60, .45)
meta.ave.cor(.05, n, cor, 0, bystudy = TRUE)
#>         Estimate      SE     LL     UL
#> Average    0.525 0.05113 0.4177 0.6179
#> Study 1    0.400 0.11431 0.1507 0.6015
#> Study 2    0.650 0.04201 0.5594 0.7252
#> Study 3    0.600 0.08000 0.4171 0.7362
#> Study 4    0.450 0.13677 0.1374 0.6811

# Should return:
#         Estimate      SE     LL     UL
# Average    0.525 0.05113 0.4177 0.6179
# Study 1    0.400 0.11431 0.1507 0.6015
# Study 2    0.650 0.04201 0.5594 0.7252
# Study 3    0.600 0.08000 0.4171 0.7362
# Study 4    0.450 0.13677 0.1374 0.6811

n <- c(150, 125, 80)
cor <- c(.396, .454, .427)
meta.ave.cor(.05, n, cor, 2, bystudy = TRUE)
#>         Estimate      SE     LL     UL
#> Average   0.4257 0.04603 0.3314 0.5115
#> Study 1   0.3960 0.06954 0.2507 0.5239
#> Study 2   0.4540 0.07187 0.3012 0.5841
#> Study 3   0.4270 0.09318 0.2259 0.5932

# Should return:
#         Estimate      SE     LL     UL
# Average   0.4257 0.04603 0.3314 0.5115

```
