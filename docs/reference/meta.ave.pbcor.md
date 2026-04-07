# Confidence interval for an average point-biserial correlation

Computes the estimate, standard error, and confidence interval for an
average point-biserial correlation from two or more studies. Two types
of point-biserial correlations can be meta-analyzed. One type uses an
unweighted variance and is appropriate in 2-group experimental designs.
The other type uses a weighted variance and is appropriate in 2-group
nonexperimental designs with simple random sampling (but not stratified
random sampling) within each study. This function requires all
point-biserial correlations to be of the same type. Use the
meta.ave.cor.gen function to meta-analyze any combination of biserial
correlation types.

For more details, see Chapter 2 of Bonett (2021, Volume 5).

## Usage

``` r
meta.ave.pbcor(alpha, m1, m2, sd1, sd2, n1, n2, type, bystudy = TRUE)
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

Bonett DG (2020). “Point-biserial correlation: Interval estimation,
hypothesis testing, meta-analysis, and sample size determination.”
*British Journal of Mathematical and Statistical Psychology*,
**73**(S1), 113–144. ISSN 0007-1102,
[doi:10.1111/bmsp.12189](https://doi.org/10.1111/bmsp.12189) .

Bonett DG (2021). *Statistical Methods for Psychologists, Vol 1-5,
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
m1 <- c(21.9, 23.1, 19.8)
m2 <- c(16.1, 17.4, 15.0)
sd1 <- c(3.82, 3.95, 3.67)
sd2 <- c(3.21, 3.30, 3.02)
n1 <- c(40, 30, 24)
n2 <- c(40, 28, 25)
meta.ave.pbcor(.05, m1, m2, sd1, sd2, n1, n2, 2, bystudy = TRUE)
#>         Estimate      SE     LL     UL
#> Average   0.6109 0.04358 0.5183 0.6893
#> Study 1   0.6350 0.06061 0.4915 0.7336
#> Study 2   0.6165 0.07458 0.4353 0.7346
#> Study 3   0.5812 0.08863 0.3648 0.7197

# Should return:
#         Estimate      SE     LL     UL
# Average   0.6109 0.04358 0.5183 0.6893
# Study 1   0.6350 0.06061 0.4915 0.7336
# Study 2   0.6165 0.07458 0.4353 0.7346
# Study 3   0.5812 0.08863 0.3648 0.7197

m1 <- c(41.2, 43.2, 49.1, 40.8)
m2 <- c(36.4, 37.1, 35.9, 31.4)
sd1 <- c(4.92, 4.75, 4.87, 5.01)
sd2 <- c(4.35, 4.24, 4.12, 4.87)
n1 <- c(42, 58, 62, 39)
n2 <- c(67, 70, 84, 45)
meta.ave.pbcor(.05, m1, m2, sd1, sd2, n1, n2, 1, bystudy = FALSE)
#>         Estimate      SE     LL     UL
#> Average   0.6358 0.02701 0.5798 0.6857
# Should return:
#         Estimate      SE     LL    UL
# Average   0.6358 0.02701 0.5798 0.6857

```
