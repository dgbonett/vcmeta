# Computes the standard error for a log odds ratio

Computes a log odds ratio and its standard error using the frequency
counts and sample sizes in a 2-group design. These frequency counts and
sample sizes can be obtained from a 2 x 2 contingency table. The log odd
ratio and its standard error are computed using a .5 addition to each
frequency count of the 2 x 2 contingency table. This function is useful
in a meta-analysis of odds ratios where some studies report the sample
odds ratio and its standard error and other studies only report the
frequency counts for a 2 x 2 contingency table. The log odds ratio and
standard error output from this function can be used as input in the
[meta.ave.gen.log](https://dgbonett.github.io/vcmeta/reference/meta.ave.gen.log.md)
function.

For more details, see Chapter 1 of Bonett (2021, Volume 5)

## Usage

``` r
se.oddsratio(f1, n1, f2, n2)
```

## Arguments

- f1:

  number of participants who have the outcome in group 1

- n1:

  sample size for group 1

- f2:

  number of participants who have the outcome in group 2

- n2:

  sample size for group 2

## Value

Returns a one-row matrix:

- Estimate - estimated log odds ratio

- SE - standard error

## References

Bonett DG, Price RM (2015). “Varying coefficient meta-analysis methods
for odds ratios and risk ratios.” *Psychological Methods*, **20**(3),
394–406. ISSN 1939-1463,
[doi:10.1037/met0000032](https://doi.org/10.1037/met0000032) .

Bonett DG (2021). *Statistical Methods for Psychologists, Vol 1-5,
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
se.oddsratio(36, 50, 21, 50)
#>                  Estimate        SE
#> Log odds ratio:  1.239501 0.4204435

# Should return: 
#                  Estimate        SE
# Log odds ratio:  1.239501 0.4204435

```
