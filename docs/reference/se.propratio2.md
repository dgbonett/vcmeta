# Computes the estimate and standard error for a 2-group log proportion ratio

Computes the Price-Bonett standard error of a 2-group proportion ratio
using the frequency count and sample size for each group. The log
proportion ratio and standard error output from this function can be
used as input in the
[meta.ave.gen.log](https://dgbonett.github.io/vcmeta/reference/meta.ave.gen.log.md)
function in applications where compatible proportion ratios from a
combination of 2-group and paired-samples studies are used in the
meta-analysis. If the proportions in each group are small (less than
.1), proportion ratios may be compatible with odds ratios and then the
[meta.ave.gen.log](https://dgbonett.github.io/vcmeta/reference/meta.ave.gen.log.md)
function could be used to meta-analyze any combination of log proportion
ratios and log odds ratios.

For more details, see Chapter 1 of Bonett (2021, Volume 5)

## Usage

``` r
se.propratio2(f1, f2, n1, n2)
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

## Value

Returns a one-row matrix:

- Estimate - estimated log proportion ratio

- SE - standard error

## References

Price RM, Bonett DG (2008). “Confidence intervals for a ratio of two
independent binomial proportions.” *Statistics in Medicine*, **27**(26),
5497–5508. ISSN 02776715,
[doi:10.1002/sim.3376](https://doi.org/10.1002/sim.3376) .

Bonett DG (2021). *Statistical Methods for Psychologists, Vol 1-5,
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
se.propratio2(31, 16, 40, 40)
#>                         Estimate        SE
#> Log proportion ratio:  0.6539265 0.2136218

# Should return:
#                         Estimate        SE
# Log proportion ratio:  0.6539265 0.2136218

```
