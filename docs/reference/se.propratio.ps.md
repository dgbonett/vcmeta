# Computes the estimate and standard error for a paired-samples log proportion ratio

Computes a large-sample standard error of a paired-samples log
proportion ratio using the frequency counts from a 2 x 2 contingency
table. The log proportion ratio and standard error output from this
function can be used as input in the
[meta.ave.gen.log](https://dgbonett.github.io/vcmeta/reference/meta.ave.gen.log.md)
function in applications where compatible proportion ratios from a
combination of 2-group and paired-samples studies are used in the
meta-analysis.

For more details, see Chapter 1 of Bonett (2021, Volume 5)

## Usage

``` r
se.propratio.ps(f00, f01, f10, f11)
```

## Arguments

- f00:

  number of participants with y = 0 and x = 0

- f01:

  number of participants with y = 0 and x = 1

- f10:

  number of participants with y = 1 and x = 0

- f11:

  number of participants with y = 1 and x = 1

## Value

Returns a one-row matrix:

- Estimate - estimated log proportion ratio

- SE - standard error

## References

Bonett DG (2021). *Statistical Methods for Psychologists, Vol 1-5,
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
se.propratio.ps(16, 64, 5, 15)
#>                         Estimate        SE
#> Log proportion ratio:  -1.373716 0.2089758

# Should return:
#                         Estimate         SE
# Log proportion ratio:  -1.373716  0.2089758

```
