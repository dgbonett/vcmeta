# Confidence interval for an average slope coefficient in a general linear model or a path model.

Computes the estimate, standard error, and confidence interval for an
average slope coefficient in a general linear model (ANOVA, ANCOVA,
multiple regression) or a path model from two or more studies. A
Satterthwaite adjustment to the degrees of freedom is used to improve
the accuracy of the confidence interval for the average slope.

For more details, see Chapter 2 of Bonett (2021, Volume 5).

## Usage

``` r
meta.ave.path(alpha, n, slope, se, s, bystudy = TRUE)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- n:

  vector of sample sizes

- slope:

  vector of slope estimates

- se:

  vector of slope standard errors

- s:

  number of predictors of the response variable

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
 
n <- c(75, 85, 250, 160)
slope <- c(1.57, 1.38, 1.08, 1.25)
se <- c(.658, .724, .307, .493)
meta.ave.path(.05, n, slope, se, 2, bystudy = TRUE)
#>         Estimate        SE          LL       UL     df
#> Average     1.32 0.2844334  0.75994525 1.880055 263.18
#> Study 1     1.57 0.6580000  0.25830097 2.881699  72.00
#> Study 2     1.38 0.7240000 -0.06026664 2.820267  82.00
#> Study 3     1.08 0.3070000  0.47532827 1.684672 247.00
#> Study 4     1.25 0.4930000  0.27623174 2.223768 157.00

#  Should return:
#          Estimate        SE          LL       UL     df
#  Average     1.32 0.2844334  0.75994528 1.880055 263.18
#  Study 1     1.57 0.6580000  0.25830097 2.881699  72.00
#  Study 2     1.38 0.7240000 -0.06026664 2.820267  82.00
#  Study 3     1.08 0.3070000  0.47532827 1.684672 247.00
#  Study 4     1.25 0.4930000  0.27623174 2.223768 157.00

```
