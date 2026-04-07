# Confidence interval for an average G-index agreement coefficient

Computes the estimate, standard error, and confidence interval for an
average G-index of agreement from two or more studies. This function
assumes that two raters each provide a dichotomous rating to a sample of
objects. As a measure of agreement, the G-index is usually preferred to
Cohen's kappa.

For more details, see Chapter 2 of Bonett (2021, Volume 5).

## Usage

``` r
meta.ave.agree(alpha, f11, f12, f21, f22, bystudy = TRUE)
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

Bonett DG (2022). “Statistical inference for G-indices of agreement.”
*Journal of Educational and Behavioral Statistics*, **47**(4), 438–458.
[doi:10.3102/10769986221088561](https://doi.org/10.3102/10769986221088561)
.

Bonett DG (2021). *Statistical Methods for Psychologists, Vol 1-5,
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
f11 <- c(43, 56, 49)
f12 <- c(7, 2, 9)
f21 <- c(3, 5, 5)
f22 <- c(37, 54, 39)
meta.ave.agree(.05, f11, f12, f21, f22, bystudy = TRUE)
#>         Estimate      SE     LL     UL
#> Average   0.7843 0.03540 0.7149 0.8537
#> Study 1   0.7447 0.06884 0.6098 0.8796
#> Study 2   0.8512 0.04771 0.7577 0.9447
#> Study 3   0.6981 0.06954 0.5618 0.8344

# Should return:
#        Estimate      SE     LL     UL
# Average  0.7843 0.03540 0.7149 0.8537
# Study 1  0.7447 0.06884 0.6098 0.8796
# Study 2  0.8512 0.04771 0.7577 0.9447
# Study 3  0.6981 0.06954 0.5618 0.8344

```
