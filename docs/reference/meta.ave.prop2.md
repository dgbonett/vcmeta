# Confidence interval for an average proportion difference in 2-group studies

Computes the estimate, standard error, and confidence interval for an
average proportion difference from two or more studies.

For more details, see Chapter 2 of Bonett (2021, Volume 5).

## Usage

``` r
meta.ave.prop2(alpha, f1, f2, n1, n2, bystudy = TRUE)
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

## References

Bonett DG, Price RM (2014). “Meta-analysis methods for risk
differences.” *British Journal of Mathematical and Statistical
Psychology*, **67**(3), 371–387. ISSN 00071102,
[doi:10.1111/bmsp.12024](https://doi.org/10.1111/bmsp.12024) .

Bonett DG (2021). *Statistical Methods for Psychologists, Vol 1-5,
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
n1 <- c(204, 201, 932, 130, 77)
n2 <- c(106, 103, 415, 132, 83)
f1 <- c(24, 40, 93, 14, 5)
f2 <- c(12, 9, 28, 3, 1)
meta.ave.prop2(.05, f1, f2, n1, n2, bystudy = TRUE)
#>             Estimate         SE            LL         UL
#> Average 0.0567907589 0.01441216  2.854345e-02 0.08503807
#> Study 1 0.0009888529 0.03870413 -7.486985e-02 0.07684756
#> Study 2 0.1067323481 0.04018243  2.797623e-02 0.18548847
#> Study 3 0.0310980338 0.01587717 -2.064379e-05 0.06221671
#> Study 4 0.0837856174 0.03129171  2.245499e-02 0.14511624
#> Study 5 0.0524199553 0.03403926 -1.429577e-02 0.11913568

# Should return:
#             Estimate         SE            LL         UL
# Average 0.0567907589 0.01441216  2.854345e-02 0.08503807
# Study 1 0.0009888529 0.03870413 -7.486985e-02 0.07684756
# Study 2 0.1067323481 0.04018243  2.797623e-02 0.18548847
# Study 3 0.0310980338 0.01587717 -2.064379e-05 0.06221671
# Study 4 0.0837856174 0.03129171  2.245499e-02 0.14511624
# Study 5 0.0524199553 0.03403926 -1.429577e-02 0.11913568

```
