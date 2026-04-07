# Confidence interval for an average Cronbach alpha reliability

Computes the estimate, standard error, and confidence interval for an
average Cronbach reliability coefficient from two or more studies.

For more details, see Chapter 2 of Bonett (2021, Volume 5).

## Usage

``` r
meta.ave.cronbach(alpha, n, rel, r, bystudy = TRUE)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- n:

  vector of sample sizes

- rel:

  vector of sample reliabilities

- r:

  number of measurements (e.g., items) used to compute each reliability

- bystudy:

  logical to also return each study estimate (TRUE) or not

## Value

Returns a matrix. The first row is the average estimate across all
studies. If bystudy is TRUE, there is 1 additional row for each study.
The matrix has the following columns:

- Estimate - estimated effect size

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
n <- c(583, 470, 546, 680)
rel <- c(.91, .89, .90, .89)
meta.ave.cronbach(.05, n, rel, 10, bystudy = TRUE)
#>         Estimate      SE     LL     UL
#> Average   0.8975 0.00326 0.8911 0.9039
#> Study 1   0.9100 0.00557 0.8986 0.9204
#> Study 2   0.8900 0.00758 0.8744 0.9041
#> Study 3   0.9000 0.00639 0.8869 0.9119
#> Study 4   0.8900 0.00630 0.8771 0.9018

# Should return:
#         Estimate      SE     LL     UL
# Average   0.8975 0.00326 0.8911 0.9039
# Study 1   0.9100 0.00557 0.8986 0.9204
# Study 2   0.8900 0.00758 0.8744 0.9041
# Study 3   0.9000 0.00639 0.8869 0.9119
# Study 4   0.8900 0.00630 0.8771 0.9018

```
