# Confidence interval for an average correlation of any type

Computes the estimate, standard error, and confidence interval for an
average correlation. Any type of correlation can be used (e.g., Pearson,
Spearman, semipartial, factor correlation, gamma coefficient, Somers d
coefficient, tetrachoric, point-biserial, biserial, correlation between
latent factors, etc.). Each study should have the same type of
correlation. If different types of correlations are used, they are
assumed to be compatible.

For more details, see Chapter 2 of Bonett (2021, Volume 5).

## Usage

``` r
meta.ave.cor.gen(alpha, cor, se, bystudy = TRUE)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- cor:

  vector of estimated correlations

- se:

  vector of standard errors

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
cor <- c(.396, .454, .409, .502, .350)
se <- c(.104, .064, .058, .107, .086)
meta.ave.cor.gen(.05, cor, se, bystudy = TRUE)
#>         Estimate      SE     LL     UL
#> Average   0.4222 0.03853 0.3439 0.4947
#> Study 1   0.3960 0.10400 0.1753 0.5788
#> Study 2   0.4540 0.06400 0.3201 0.5701
#> Study 3   0.4090 0.05800 0.2894 0.5160
#> Study 4   0.5020 0.10700 0.2651 0.6817
#> Study 5   0.3500 0.08600 0.1716 0.5061

# Should return:
#         Estimate      SE     LL     UL
# Average   0.4222 0.03853 0.3439 0.4947
# Study 1   0.3960 0.10400 0.1753 0.5788
# Study 2   0.4540 0.06400 0.3201 0.5701
# Study 3   0.4090 0.05800 0.2894 0.5160
# Study 4   0.5020 0.10700 0.2651 0.6817
# Study 5   0.3500 0.08600 0.1716 0.5061

```
