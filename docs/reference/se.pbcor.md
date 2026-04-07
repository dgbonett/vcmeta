# Computes the standard error for a point-biserial correlation

Computes a point-biserial correlation and its standard error for two
types of point-biserial correlations in 2-group designs using the
estimated means, estimated standard deviations, and sample sizes.
Equality of variances is not assumed. One type of point-biserial
correlation uses an unweighted average of variances and is recommended
for 2-group experimental designs. The other type of point-biserial
correlation uses a weighted average of variances and is recommended for
2-group nonexperimental designs with simple random sampling (but not
stratified random sampling). This function is useful in a meta-analysis
of compatible point-biserial correlations where some studies used a
2-group experimental design and other studies used a 2-group
nonexperimental design. The effect size estimate and standard error
output from this function can be used as input in the
[meta.ave.cor.gen](https://dgbonett.github.io/vcmeta/reference/meta.ave.cor.gen.md)
function.

For more details, see Chapter 1 of Bonett (2021, Volume 5)

## Usage

``` r
se.pbcor(m1, m2, sd1, sd2, n1, n2, type)
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

- type:

  - set to 1 for weighted variance average

  - set to 2 for unweighted variance average

## Value

Returns a one-row matrix:

- Estimate - estimated point-biserial correlation

- SE - standard error

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
se.pbcor(21.9, 16.1, 3.82, 3.21, 40, 40, 1)
#>                              Estimate      SE
#> Point-biserial correlation:     0.635 0.05981

#  Should return: 
#                              Estimate      SE
#  Point-biserial correlation:    0.635 0.05981

```
