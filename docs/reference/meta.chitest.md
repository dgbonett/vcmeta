# Computes a chi-square test of effect-size homogeneity

Computes a chi-square test of effect size homogeneity and p-value using
effect-size estimates and their standard errors from two or more
studies. This test should not be used to justify the use of a constant
coeffient (fixed-effect) meta-analysis.

## Usage

``` r
meta.chitest(est, se)
```

## Arguments

- est:

  vector of effect-size estimates

- se:

  vector of effect-size standard errors

## Value

Returns a one-row matrix:

- Q - chi-square test statitic

- df - degrees of freedom

- p - p-value

## References

Borenstein M, Hedges LV, Higgins JP, Rothstein HR (2009). *Introduction
to meta-analysis*. Wiley, New York.

## Examples

``` r
est <- c(.297, .324, .281, .149) 
se <- c(.082, .051, .047, .094)
meta.chitest(est, se)
#>      Q df     p
#>  2.707  3 0.439

# Should return:
#      Q df     p
#  2.707  3 0.439

```
