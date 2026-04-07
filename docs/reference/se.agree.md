# Computes the estimate and standard error for a G-index of agreement

Computes an adjusted G-index of agreement for two raters using the
number of objects rated in agreement, the sample size (number of
objects), and the number of rating categories. The adjustment for the
G-index and its standard error is optimized for the the number of
planned studies in the meta-analysis. The G-index and standard error
output from this function can be used as input in the
[meta.ave.gen](https://dgbonett.github.io/vcmeta/reference/meta.ave.gen.md),
[meta.lc.gen](https://dgbonett.github.io/vcmeta/reference/meta.lc.gen.md),
and
[meta.lm.gen](https://dgbonett.github.io/vcmeta/reference/meta.lm.gen.md)
functions. The G-index is usually preferred to Cohen's kappa (see
Bonett, 2022).

For more details, see Chapter 1 of Bonett (2021, Volume 5)

## Usage

``` r
se.agree(f, n, k, m)
```

## Arguments

- f:

  number of objects rated in agreement

- n:

  sample size (number of objects)

- k:

  number of rating categories

- m:

  number of studies in planned meta-analysis

## Value

Returns a one-row matrix:

- MLE - maximum likelihood estimate of G-index

- Estimate - adjusted estimate of G-index for meta-analysis

- SE - standard error of adjusted estimate for meta-analysis

## References

Bonett DG (2022). “Statistical inference for G-indices of agreement.”
*Journal of Educational and Behavioral Statistics*, **47**(4), 438–458.
[doi:10.3102/10769986221088561](https://doi.org/10.3102/10769986221088561)
.

Bonett DG (2021). *Statistical Methods for Psychologists, Vol 1-5,
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
se.agree(42, 50, 3, 4)
#>            MLE Estimate      SE
#> G-index:  0.76     0.75 0.06391

# Should return: 
#             MLE  Estimate      SE
# G-index:   0.76      0.75 0.06391 

```
