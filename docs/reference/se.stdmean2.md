# Computes the standard error for a 2-group standardized mean difference

Computes the standard error of a 2-group standardized mean difference
using the sample sizes and the estimated means and standard deviations.
Use the square root average variance standardizer (stdzr = 0) for
2-group experimental designs. Use the square root weighted variance
standardizer (stdzr = 3) for 2-group nonexperimental designs with simple
random sampling. The single-group standardizers (stdzr = 1 and stdzr =
2) can be used with either 2-group experimental or nonexperimental
designs. The effect size estimate and standard error output from this
function can be used as input in the
[meta.ave.gen](https://dgbonett.github.io/vcmeta/reference/meta.ave.gen.md),
[meta.lc.gen](https://dgbonett.github.io/vcmeta/reference/meta.lc.gen.md),
and
[meta.lm.gen](https://dgbonett.github.io/vcmeta/reference/meta.lm.gen.md)
functions in applications where compatible standardized mean differences
from a combination of 2-group and paired-samples experiments are used in
the meta-analysis. Equality of variances is not assumed.

For more details, see Chapter 1 of Bonett (2021, Volume 5)

## Usage

``` r
se.stdmean2(m1, m2, sd1, sd2, n1, n2, stdzr)
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

- stdzr:

  - set to 0 for square root average variance standardizer

  - set to 1 for group 1 SD standardizer

  - set to 2 for group 2 SD standardizer

  - set to 3 for square root weighted variance standardizer

## Value

Returns a one-row matrix:

- Estimate - estimated standardized mean difference

- SE - standard error

## References

Bonett DG (2009). “Meta-analytic interval estimation for standardized
and unstandardized mean differences.” *Psychological Methods*,
**14**(3), 225–238. ISSN 1939-1463,
[doi:10.1037/a0016619](https://doi.org/10.1037/a0016619) .

Bonett DG (2021). *Statistical Methods for Psychologists, Vol 1-5,
https://dgbonett.sites.ucsc.edu/*.

## See also

[se.cohen](https://dgbonett.github.io/vcmeta/reference/se.cohen.md)

## Examples

``` r
se.stdmean2(21.9, 16.1, 3.82, 3.21, 40, 40, 0)
#>                                Estimate     SE
#> Standardized mean difference:    1.6439 0.2629

# Should return: 
#                               Estimate      SE
# Standardized mean difference:   1.6439 0.26290

se.stdmean2(21.9, 16.1, 3.82, 3.21, 31, 49, 3)
#>                                Estimate      SE
#> Standardized mean difference:    1.6776 0.27573

# Should return: 
#                               Estimate      SE
# Standardized mean difference:   1.6776 0.27573

```
