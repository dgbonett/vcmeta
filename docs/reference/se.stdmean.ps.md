# Computes the standard error for a paired-samples standardized mean difference

Computes the standard error of a paired-samples standardized mean
difference using the sample size, estimated means, estimated standard
deviations, and estimated correlation. The effect size estimate and
standard error output from this function can be used as input in the
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
se.stdmean.ps(m1, m2, sd1, sd2, cor, n, stdzr)
```

## Arguments

- m1:

  estimated mean for measurement 1

- m2:

  estimated mean for measurement 2

- sd1:

  estimated standard deviation for measurement 1

- sd2:

  estimated standard deviation for measurement 2

- cor:

  estimated correlation for measurements 1 and 2

- n:

  sample size

- stdzr:

  - set to 0 for square root average variance standardizer

  - set to 1 for measurement 1 SD standardizer

  - set to 2 for measurement 2 SD standardizer

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

## Examples

``` r
se.stdmean.ps(23.9, 25.1, 1.76, 2.01, .78, 25, 0)
#>                                Estimate      SE
#> Standardized mean difference:   -0.6352 0.16029

# Should return: 
#                                 Estimate      SE
# Standardized mean difference:   -0.6352 0.16029

se.stdmean.ps(23.9, 25.1, 1.76, 2.01, .78, 25, 1)
#>                                Estimate      SE
#> Standardized mean difference:   -0.6818 0.17738

# Should return:
#                                Estimate      SE
# Standardized mean difference:   -0.6818 0.17738

```
