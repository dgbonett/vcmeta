# Computes the standard error for a biserial-phi correlation

Computes an estimate of a biserial-phi correlation and its standard
error using the frequency counts from a 2 x 2 contingency table where
one variable is naturally dichotomous and the other variable is
artifically dichotomous. The biserial-phi correlation approximates a
point-biserial correlation between the naturally dichotomous variable
and the unobserved quantitative variable that was measured on a
dichotomous scale. A biserial-phi correlation could be compatible with a
point-biserial correlation in a meta-analysis. The biserial-phi estimate
and the standard error from this function can be used as input in the
[meta.ave.cor.gen](https://dgbonett.github.io/vcmeta/reference/meta.ave.cor.gen.md)
function in a meta-analysis where a point-biserial correlation has been
obtained in some studies and a biserial-phi correlation has been
obtained in other studies.

For more details, see Chapter 1 of Bonett (2021, Volume 5)

## Usage

``` r
se.biphi(f1, f2, n1, n2)
```

## Arguments

- f1:

  number of participants in group 1 who have the attribute

- f2:

  number of participants in group 2 who have the attribute

- n1:

  sample size for group 1

- n2:

  sample size for group 2

## Value

Returns a 1-row matrix. The columns are:

- Estimate - estimated biserial-phi correlation

- SE - standard error

## References

Bonett DG (2021). *Statistical Methods for Psychologists, Vol 1-5,
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
se.biphi(34, 22, 50, 50)
#>                Estimate      SE
#> Biserial-phi:    0.2754 0.10746

# Should return:
#               Estimate      SE 
# Biserial-phi:   0.2754 0.10746

```
