# Confidence interval for a subgroup difference in average Cronbach reliabilities

Computes the estimate, standard error, and confidence interval for a
difference in average Cronbach reliability coefficients for two mutually
exclusive subgroups of studies. Each set can have one or more studies.
The number of measurements used to compute the sample reliablity
coefficient is assumed to be the same for all studies.

For more details, see Section 3.3 of Bonett (2021, Volume 5).

## Usage

``` r
meta.sub.cronbach(alpha, n, rel, r, group)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- n:

  vector of sample sizes

- rel:

  vector of estimated Cronbach reliabilities

- r:

  number of measurements (e.g., items)

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

Bonett DG (2010). “Varying coefficient meta-analytic methods for alpha
reliability.” *Psychological Methods*, **15**(4), 368–385. ISSN
1939-1463, [doi:10.1037/a0020142](https://doi.org/10.1037/a0020142) .

Bonett DG, Wright TA (2015). “Cronbach's alpha reliability: Interval
estimation, hypothesis testing, and sample size planning.” *Journal of
Organizational Behavior*, **36**(1), 3–15. ISSN 08943796,
[doi:10.1002/job.1960](https://doi.org/10.1002/job.1960) .

Bonett DG (2021). *Statistical Methods for Psychologists, Vol 1-5,
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
n <- c(120, 170, 150, 135)
rel <- c(.891, .873, .734, .715)
group <- c(1, 1, 2, 2)
meta.sub.cronbach(.05, n, rel, 10, group)
#>                Estimate      SE     LL     UL
#> Set A:           0.8820 0.01052 0.8605 0.9016
#> Set B:           0.7245 0.02474 0.6738 0.7706
#> Set A - Set B:   0.1575 0.02689 0.1066 0.2119

# Should return: 
#                Estimate      SE     LL     UL
# Set A:           0.8820 0.01052 0.8605 0.9016
# Set B:           0.7245 0.02474 0.6738 0.7706
# Set A - Set B:   0.1575 0.02689 0.1066 0.2119

```
