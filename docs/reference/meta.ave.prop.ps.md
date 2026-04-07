# Confidence interval for an average proportion difference in paired-samples studies

Computes the estimate, standard error, and confidence interval for an
average proportion difference from two or more studies.

For more details, see Chapter 2 of Bonett (2021, Volume 5).

## Usage

``` r
meta.ave.prop.ps(alpha, f11, f12, f21, f22, bystudy = TRUE)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- f11:

  vector of frequency counts in cell 1,1

- f12:

  vector of frequency counts in cell 1,2

- f21:

  vector of frequency counts in cell 2,1

- f22:

  vector of frequency counts in cell 2,2

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

Bonett DG, Price RM (2012). “Adjusted wald confidence interval for a
difference of binomial proportions based on paired data.” *Journal of
Educational and Behavioral Statistics*, **37**(4), 479–488. ISSN
1076-9986,
[doi:10.3102/1076998611411915](https://doi.org/10.3102/1076998611411915)
.

Bonett DG (2021). *Statistical Methods for Psychologists, Vol 1-5,
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
f11 <- c(17, 28, 19)
f12 <- c(43, 56, 49)
f21 <- c(3, 5, 5)
f22 <- c(37, 54, 39)
meta.ave.prop.ps(.05, f11, f12, f21, f22, bystudy = TRUE)
#>          Estimate         SE        LL        UL
#> Average 0.3809573 0.03000016 0.3221581 0.4397565
#> Study 1 0.3921569 0.05573055 0.2829270 0.5013867
#> Study 2 0.3517241 0.04629537 0.2609869 0.4424614
#> Study 3 0.3859649 0.05479300 0.2785726 0.4933572

# Should return:
#          Estimate         SE        LL        UL
# Average 0.3809573 0.03000016 0.3221581 0.4397565
# Study 1 0.3921569 0.05573055 0.2829270 0.5013867
# Study 2 0.3517241 0.04629537 0.2609869 0.4424614
# Study 3 0.3859649 0.05479300 0.2785726 0.4933572

```
