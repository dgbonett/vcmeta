# Computes the standard error for the average of 2-group mean differences from two parallel measurement response variables in the same sample

In a study that reports a 2-group mean difference for two response
variables that satisfy the conditions of parallel measurments, this
function can be used to compute the standard error of the average of the
two mean differences using the two estimated means, estimated standard
deviations, estimated within-group correlation between the two response
variables, and the two sample sizes. The average mean difference and
standard error output from this function can then be used as input in
the
[meta.ave.gen](https://dgbonett.github.io/vcmeta/reference/meta.ave.gen.md),
[meta.lc.gen](https://dgbonett.github.io/vcmeta/reference/meta.lc.gen.md),
and
[meta.lm.gen](https://dgbonett.github.io/vcmeta/reference/meta.lm.gen.md)
functions in a meta-analysis where some studies have used one of the two
parallel response variables and other studies have used the other
parallel response variable. Equality of variances is not assumed.

For more details, see Chapter 1 of Bonett (2021, Volume 5)

## Usage

``` r
se.ave.mean2.dep(m1A, m2A, sd1A, sd2A, m1B, m2B, sd1B, sd2B, rAB, n1, n2)
```

## Arguments

- m1A:

  estimated mean for variable A in group 1

- m2A:

  estimated mean for variable A in group 2

- sd1A:

  estimated standard deviation for variable A in group 1

- sd2A:

  estimated standard deviation for variable A in group 2

- m1B:

  estimated mean for variable B in group 1

- m2B:

  estimated mean for variable B in group 2

- sd1B:

  estimated standard deviation for variable B in group 1

- sd2B:

  estimated standard deviation for variable B in group 2

- rAB:

  estimated within-group correlation between variables A and B

- n1:

  sample size for group 1

- n2:

  sample size for group 2

## Value

Returns a one-row matrix:

- Estimate - estimated average mean difference

- SE - standard error

- VAR(A) - variance of mean difference for variable A

- VAR(B) - variance of mean difference for variable B

- COV(A,B) - covariance of mean differences for variables A and B

## References

Bonett DG (2021). *Statistical Methods for Psychologists, Vol 1-5,
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
se.ave.mean2.dep(21.9, 16.1, 3.82, 3.21, 24.8, 17.1, 3.57, 3.64, .785, 40, 40)
#>                           Estimate        SE    VAR(A)    VAR(B)  COV(A,B)
#> Average mean difference:      6.75 0.7526878 0.6224125 0.6498625 0.4969403

# Should return:
#                          Estimate        SE    VAR(A)    VAR(B)  COV(A,B)
# Average mean difference:     6.75 0.7526878 0.6224125 0.6498625 0.4969403

```
