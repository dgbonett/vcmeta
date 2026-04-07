# Confidence interval for an average of any parameter

Computes the estimate, standard error, and confidence interval for an
average of any type of parameter from two or more studies. Each study
should have the same type of parameter.

For more details, see Chapter 2 of Bonett (2021, Volume 5).

## Usage

``` r
meta.ave.gen(alpha, est, se, bystudy = TRUE)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- est:

  vector of parameter estimates

- se:

  vector of standard errors

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
est <- c(.022, .751, .421, .287, .052, .146, .562, .904)
se <- c(.124, .464, .102, .592, .864, .241, .252, .318)
meta.ave.gen(.05, est, se, bystudy = TRUE)
#>         Estimate        SE          LL        UL
#> Average 0.393125 0.1561622  0.08705266 0.6991973
#> Study 1 0.022000 0.1240000 -0.22103553 0.2650355
#> Study 2 0.751000 0.4640000 -0.15842329 1.6604233
#> Study 3 0.421000 0.1020000  0.22108367 0.6209163
#> Study 4 0.287000 0.5920000 -0.87329868 1.4472987
#> Study 5 0.052000 0.8640000 -1.64140888 1.7454089
#> Study 6 0.146000 0.2410000 -0.32635132 0.6183513
#> Study 7 0.562000 0.2520000  0.06808908 1.0559109
#> Study 8 0.904000 0.3180000  0.28073145 1.5272685

# Should return:
#          Estimate        SE          LL        UL
# Average  0.393125 0.1561622  0.08705266 0.6991973
# Study 1  0.022000 0.1240000 -0.22103553 0.2650355
# Study 2  0.751000 0.4640000 -0.15842329 1.6604233
# Study 3  0.421000 0.1020000  0.22108367 0.6209163
# Study 4  0.287000 0.5920000 -0.87329868 1.4472987
# Study 5  0.052000 0.8640000 -1.64140888 1.7454089
# Study 6  0.146000 0.2410000 -0.32635132 0.6183513
# Study 7  0.562000 0.2520000  0.06808908 1.0559109
# Study 8  0.904000 0.3180000  0.28073145 1.5272685

```
