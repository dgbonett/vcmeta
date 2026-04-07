# Computes the estimate and standard error for a 2-group proportion difference

This function computes the Price-Bonett standard error of a 2-group
proportion difference using the frequency counts, sample sizes, and
planned number of studies in the meta-analysis. The effect size estimate
and standard error output from this function can be used as input in the
[meta.ave.gen](https://dgbonett.github.io/vcmeta/reference/meta.ave.gen.md),
[meta.lc.gen](https://dgbonett.github.io/vcmeta/reference/meta.lc.gen.md),
and
[meta.lm.gen](https://dgbonett.github.io/vcmeta/reference/meta.lm.gen.md)
functions in applications where compatible proportion differences from a
combination of 2-group and paired-samples studies are used in the
meta-analysis.

For more details, see Chapter 1 of Bonett (2021, Volume 5)

## Usage

``` r
se.prop2(f1, f2, n1, n2, m)
```

## Arguments

- f1:

  number of participants in group 1 who have the outcome

- f2:

  number of participants in group 2 who have the outcome

- n1:

  sample size for group 1

- n2:

  sample size for group 2

- m:

  number of studies in planned meta-analysis

## Value

Returns a one-row matrix:

- Estimate - adjusted estimate of proportion difference for
  meta-analysis

- SE - standard error of adjusted estimate for meta-analysis

## References

Price RM, Bonett DG (2004). “An improved confidence interval for a
linear function of binomial proportions.” *Computational Statistics and
Data Analysis*, **45**(3), 449–456. ISSN 01679473,
[doi:10.1016/S0167-9473(03)00007-0](https://doi.org/10.1016/S0167-9473%2803%2900007-0)
.

Bonett DG, Price RM (2014). “Meta-analysis methods for risk
differences.” *British Journal of Mathematical and Statistical
Psychology*, **67**(3), 371–387. ISSN 00071102,
[doi:10.1111/bmsp.12024](https://doi.org/10.1111/bmsp.12024) .

Bonett DG (2021). *Statistical Methods for Psychologists, Vol 1-5,
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
se.prop2(31, 16, 40, 40, 5)
#>                          Estimate        SE
#> Proportion difference:  0.3712871 0.1014819

# Should return:
#                          Estimate        SE
# Proportion difference:  0.3676471 0.1011817

```
