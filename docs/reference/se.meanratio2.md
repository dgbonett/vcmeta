# Computes the standard error for a 2-group log mean ratio

Computes the standard error of a 2-group log mean ratio using the
estimated means, estimated standard deviations, and sample sizes. The
log mean estimate and standard error output from this function can be
used as input in the
[meta.ave.gen](https://dgbonett.github.io/vcmeta/reference/meta.ave.gen.md),
[meta.lc.gen](https://dgbonett.github.io/vcmeta/reference/meta.lc.gen.md),
and
[meta.lm.gen](https://dgbonett.github.io/vcmeta/reference/meta.lm.gen.md)
functions in application where compatible mean ratios from a combination
of 2-group and paired-samples experiments are used in the meta-analysis.
Equality of variances is not assumed.

For more details, see Chapter 1 of Bonett (2021, Volume 5)

## Usage

``` r
se.meanratio2(m1, m2, sd1, sd2, n1, n2)
```

## Arguments

- m1:

  estimated mean for group 1

- m2:

  estimated mean for group 2

- sd1:

  estimated standard deviation for group 1

- sd2:

  estimated standard deviation for group 2

- n1:

  sample size for group 1

- n2:

  sample size for group 2

## Value

Returns a one-row matrix:

- Estimate - estimated log mean ratio

- SE - standard error

## References

Bonett DG, Price RM (2020). “Confidence intervals for ratios of means
and medians.” *Journal of Educational and Behavioral Statistics*,
**45**(6), 750–770. ISSN 1076-9986,
[doi:10.3102/1076998620934125](https://doi.org/10.3102/1076998620934125)
.

Bonett DG (2021). *Statistical Methods for Psychologists, Vol 1-5,
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
se.meanratio2(21.9, 16.1, 3.82, 3.21, 40, 40)
#>                   Estimate       SE
#> Log mean ratio:  0.3076674 0.041886

# Should return:
#                   Estimate       SE
# Log mean ratio:  0.3076674 0.041886

```
