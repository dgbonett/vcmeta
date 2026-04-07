# Computes a Pearson correlation between paired measurements from a paired-samples t statistic

Computes the Pearson correlation between paired measurements using a
reported paired-samples t statistic and other sample information. This
correlation estimate is needed in several functions that analyze mean
differences and standardized mean differences in paired-samples studies.

## Usage

``` r
cor.from.t(m1, m2, sd1, sd2, t, n)
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

- t:

  value of paired-samples t-test

- n:

  sample size

## Value

Returns the sample Pearson correlation between the two paired
measurements

## Examples

``` r
cor.from.t(9.4, 9.8, 1.26, 1.40, 2.27, 30)
#>               Estimate
#> Correlation:    0.7415

# Should return:
#              Estimate
# Correlation:   0.7415

```
