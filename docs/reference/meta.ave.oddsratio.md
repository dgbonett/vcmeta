# Confidence interval for average odds ratio from 2-group studies

Computes the estimate, standard error, and confidence interval for a
geometric average odds ratio from two or more studies.

For more details, see Chapter 2 of Bonett (2021, Volume 5).

## Usage

``` r
meta.ave.oddsratio(alpha, f1, f2, n1, n2, bystudy = TRUE)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- f1:

  vector of group 1 frequency counts

- f2:

  vector of group 2 frequency counts

- n1:

  vector of group 1 sample sizes

- n2:

  vector of group 2 sample sizes

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

- exp(Estimate) - the exponentiated estimate

- exp(LL) - lower limit of the exponentiated confidence interval

- exp(UL) - upper limit of the exponentiated confidence interval

## References

Bonett DG, Price RM (2015). “Varying coefficient meta-analysis methods
for odds ratios and risk ratios.” *Psychological Methods*, **20**(3),
394–406. ISSN 1939-1463,
[doi:10.1037/met0000032](https://doi.org/10.1037/met0000032) .

Bonett DG (2021). *Statistical Methods for Psychologists, Vol 1-5,
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
n1 <- c(204, 201, 932, 130, 77)
n2 <- c(106, 103, 415, 132, 83)
f1 <- c(24, 40, 93, 14, 5)
f2 <- c(12, 9, 28, 3, 1)
meta.ave.oddsratio(.05, f1, f2, n1, n2, bystudy = TRUE)
#>           Estimate        SE          LL        UL exp(Estimate)   exp(LL)
#> Average 0.86211102 0.2512852  0.36960107 1.3546210      2.368155 1.4471572
#> Study 1 0.02581353 0.3700520 -0.69947512 0.7511022      1.026150 0.4968460
#> Study 2 0.91410487 0.3830515  0.16333766 1.6648721      2.494541 1.1774342
#> Study 3 0.41496672 0.2226089 -0.02133877 0.8512722      1.514320 0.9788873
#> Study 4 1.52717529 0.6090858  0.33338907 2.7209615      4.605150 1.3956902
#> Study 5 1.42849472 0.9350931 -0.40425414 3.2612436      4.172414 0.6674745
#>           exp(UL)
#> Average  3.875292
#> Study 1  2.119335
#> Study 2  5.284997
#> Study 3  2.342625
#> Study 4 15.194925
#> Study 5 26.081952

# Should return:
#           Estimate        SE          LL        UL 
# Average 0.86211102 0.2512852  0.36960107 1.3546210
# Study 1 0.02581353 0.3700520 -0.69947512 0.7511022
# Study 2 0.91410487 0.3830515  0.16333766 1.6648721
# Study 3 0.41496672 0.2226089 -0.02133877 0.8512722
# Study 4 1.52717529 0.6090858  0.33338907 2.7209615
# Study 5 1.42849472 0.9350931 -0.40425414 3.2612436
#         exp(Estimate)   exp(LL)   exp(UL)
# Average      2.368155 1.4471572  3.875292
# Study 1      1.026150 0.4968460  2.119335
# Study 2      2.494541 1.1774342  5.284997
# Study 3      1.514320 0.9788873  2.342625
# Study 4      4.605150 1.3956902 15.194925
# Study 5      4.172414 0.6674745 26.081952

```
