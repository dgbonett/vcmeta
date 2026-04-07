# Confidence interval for an average mean ratio from paired-samples studies

Computes the estimate, standard error, and confidence interval for a
geometric average mean ratio from two or more paired-samples studies. A
Satterthwaite adjustment to the degrees of freedom is used to improve
the accuracy of the confidence interval for the average effect size.
Equality of variances within or across studies is not assumed.

For more details, see Chapter 2 of Bonett (2021, Volume 5).

## Usage

``` r
meta.ave.meanratio.ps(alpha, m1, m2, sd1, sd2, cor, n, bystudy = TRUE)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- m1:

  vector of estimated means for measurement 1

- m2:

  vector of estimated means for measurement 2

- sd1:

  vector of estimated SDs for measurement 1

- sd2:

  vector of estimated SDs for measurement 2

- cor:

  vector of estimated correlations for paired measurements

- n:

  vector of sample sizes

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

- exp(Estimate) - exponentiated estimate

- exp(LL) - lower limit of the exponentiated confidence interval

- exp(UL) - upper limit of the exponentiated confidence interval

- df - degrees of freedom

## References

Bonett DG, Price RM (2020). “Confidence intervals for ratios of means
and medians.” *Journal of Educational and Behavioral Statistics*,
**45**(6), 750–770. ISSN 1076-9986,
[doi:10.3102/1076998620934125](https://doi.org/10.3102/1076998620934125)
.

Bonett DG (2021). *Statistical Methods for Psychologists, Vol 1-5,
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
m1 <- c(53, 60, 53, 57)
m2 <- c(55, 62, 58, 61)
sd1 <- c(4.1, 4.2, 4.5, 4.0)
sd2 <- c(4.2, 4.7, 4.9, 4.8)
cor <- c(.7, .7, .8, .85)
n <- c(30, 50, 30, 70)
meta.ave.meanratio.ps(.05, m1, m2, sd1, sd2, cor, n, bystudy = TRUE)
#>            Estimate          SE          LL          UL exp(Estimate)   exp(LL)
#> Average -0.05695120 0.004350863 -0.06558008 -0.04832232     0.9446402 0.9365240
#> Study 1 -0.03704127 0.010871086 -0.05927514 -0.01480740     0.9636364 0.9424474
#> Study 2 -0.03278982 0.008021952 -0.04891054 -0.01666911     0.9677419 0.9522663
#> Study 3 -0.09015110 0.009779919 -0.11015328 -0.07014892     0.9137931 0.8956968
#> Study 4 -0.06782260 0.004970015 -0.07773750 -0.05790769     0.9344262 0.9252073
#>           exp(UL)     df
#> Average 0.9528266 103.03
#> Study 1 0.9853017  29.00
#> Study 2 0.9834691  49.00
#> Study 3 0.9322550  29.00
#> Study 4 0.9437371  69.00

# Should return:
#            Estimate          SE          LL          UL
# Average -0.05695120 0.004350863 -0.06558008 -0.04832231
# Study 1 -0.03704127 0.010871086 -0.05927514 -0.01480740
# Study 2 -0.03278982 0.008021952 -0.04891054 -0.01666911
# Study 3 -0.09015110 0.009779919 -0.11015328 -0.07014892
# Study 4 -0.06782260 0.004970015 -0.07773750 -0.05790769
#         exp(Estimate)   exp(LL)   exp(UL)     df
# Average     0.9446402 0.9365240 0.9528266 103.03
# Study 1     0.9636364 0.9424474 0.9853017  29.00
# Study 2     0.9677419 0.9522663 0.9834691  49.00
# Study 3     0.9137931 0.8956968 0.9322550  29.00
# Study 4     0.9344262 0.9252073 0.9437371  69.00

```
