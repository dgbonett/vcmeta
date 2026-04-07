# Computes the standard error for a biserial correlation

Computes a biserial correlation and its standard error. A biserial
correlation can be used when one variable is quantitative and the other
variable has been artifically dichotmized. The biserial correlation
estimates the correlation between an observable quantitative variable
and an unobserved quantitative variable that is measured on a
dichotomous scale. This function requires the estimated mean, estimated
standard deviation, and samples size from each level of the dichotomized
variable. This function is useful in a meta-analysis of Pearson
correlations where some studies report a Pearson correlation and other
studies report the information needed to compute a biserial correlation.
The biserial correlation and standard error output from this function
can be used as input in the
[meta.ave.cor.gen](https://dgbonett.github.io/vcmeta/reference/meta.ave.cor.gen.md)
function.

For more details, see Chapter 1 of Bonett (2021, Volume 5)

## Usage

``` r
se.bscor(m1, m2, sd1, sd2, n1, n2)
```

## Arguments

- m1:

  estimated mean for level 1

- m2:

  estimated mean for level 2

- sd1:

  estimated standard deviation for level 1

- sd2:

  estimated standard deviation for level 2

- n1:

  sample size for level 1

- n2:

  sample size for level 2

## Value

Returns a one-row matrix:

- Estimate - estimated biserial correlation

- SE - standard error

## Details

This function computes a point-biserial correlation and its standard
error as a function of a standardized mean difference with a weighted
variance standardizer. Then the point-biserial estimate is transformed
into a biserial correlation using the traditional adjustment. The
adjustment is also applied to the point-biserial standard error to
obtain the standard error for the biserial correlation.

The biserial correlation assumes that the observed quantitative variable
and the unobserved quantitative variable have a bivariate normal
distribution. Bivariate normality is a crucial assumption underlying the
transformation of a point-biserial correlation to a biserial
correlation. Bivariate normality also implies equal variances of the
observed quantitative variable at each level of the dichotomized
variable, and this assumption is made in the computation of the standard
error.

## References

Bonett DG (2020). “Point-biserial correlation: Interval estimation,
hypothesis testing, meta-analysis, and sample size determination.”
*British Journal of Mathematical and Statistical Psychology*,
**73**(S1), 113–144. ISSN 0007-1102,
[doi:10.1111/bmsp.12189](https://doi.org/10.1111/bmsp.12189) .

Bonett DG (2021). *Statistical Methods for Psychologists, Vol 1-5,
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
se.bscor(21.9, 16.1, 3.82, 3.21, 40, 40)
#>                        Estimate      SE
#> Biserial correlation:    0.8018 0.07452

#  Should return: 
#                        Estimate      SE
#  Biserial correlation:   0.8018 0.07452

```
