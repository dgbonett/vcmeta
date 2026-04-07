# Confidence interval for a subgroup difference in average point-biserial correlations

Computes the estimate, standard error, and confidence interval for a
difference in average point-biserial correlations for two mutually
exclusive subgroups of studies. Each subgroup can have one or more
studies. Two types of point-biserial correlations can be analyzed. One
type uses an unweighted variance and is recommended for 2-group
experimental designs. The other type uses a weighted variance and is
recommended for 2-group nonexperimental designs with simple random
sampling (but not stratified random sampling) within each study.
Equality of variances within or across studies is not assumed.

For more details, see Section 3.3 of Bonett (2021, Volume 5).

## Usage

``` r
meta.sub.pbcor(alpha, m1, m2, sd1, sd2, n1, n2, type, group)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- m1:

  vector of estimated means for group 1

- m2:

  vector of estimated means for group 2

- sd1:

  vector of estimated SDs for group 1

- sd2:

  vector of estimated SDs for group 2

- n1:

  vector of group 1 sample sizes

- n2:

  vector of group 2 sample sizes

- type:

  - set to 1 for weighted variance

  - set to 2 for unweighted variance

- group:

  vector of group indicators:

  - 1 for set A

  - 2 for set B

  - 0 to ignore

## Value

Returns a matrix with three rows:

- Row 1 - estimate for Set A

- Row 2 - estimate for Set B

- Row 3 - estimate for difference, Set A - Set B

The columns are:

- Estimate - estimated average correlation or difference

- SE - standard error

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

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
m1 <- c(45.1, 39.2, 36.3, 34.5)
m2 <- c(30.0, 35.1, 35.3, 36.2)
sd1 <- c(10.7, 10.5, 9.4, 11.5)
sd2 <- c(12.3, 12.0, 10.4, 9.6)
n1 <- c(40, 20, 50, 25)
n2 <- c(40, 20, 48, 26)
group <- c(1, 1, 2, 2)
meta.sub.pbcor(.05,  m1, m2, sd1, sd2, n1, n2, 2, group)
#>                Estimate      SE      LL     UL
#> Set A:           0.3634 0.08553  0.1855 0.5182
#> Set B:          -0.0148 0.08741 -0.1840 0.1553
#> Set A - Set B:   0.3782 0.12229  0.1321 0.6076

# Should return:
#                Estimate      SE      LL     UL
# Set A:           0.3634 0.08553  0.1855 0.5182
# Set B:          -0.0148 0.08741 -0.1840 0.1553
# Set A - Set B:   0.3782 0.12229  0.1321 0.6076

```
