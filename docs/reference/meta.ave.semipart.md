# Confidence interval for an average semipartial correlation

Computes the estimate, standard error, and confidence interval for an
average semipartial correlation from two or more studies.

For more details, see Chapter 2 of Bonett (2021, Volume 5).

## Usage

``` r
meta.ave.semipart(alpha, n, cor, r2, bystudy = TRUE)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- n:

  vector of sample sizes

- cor:

  vector of estimated semipartial correlations

- r2:

  vector of squared multiple correlations for a model that includes the
  IV and all control variables

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

Bonett DG (2021). *Statistical Methods for Psychologists, Vol 1-5,
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
n <- c(128, 97, 210, 217)
cor <- c(.35, .41, .44, .39)
r2 <- c(.29, .33, .36, .39)
meta.ave.semipart(.05, n, cor, r2, bystudy = TRUE)
#>         Estimate      SE     LL     UL
#> Average   0.3975 0.03221 0.3326 0.4587
#> Study 1   0.3500 0.07175 0.2023 0.4821
#> Study 2   0.4100 0.07886 0.2447 0.5521
#> Study 3   0.4400 0.05147 0.3338 0.5351
#> Study 4   0.3900 0.05085 0.2860 0.4849

# Should return:
#         Estimate      SE     LL     UL
# Average   0.3975 0.03221 0.3326 0.4587
# Study 1   0.3500 0.07175 0.2023 0.4821
# Study 2   0.4100 0.07886 0.2447 0.5521
# Study 3   0.4400 0.05147 0.3338 0.5351
# Study 4   0.3900 0.05085 0.2860 0.4849

```
