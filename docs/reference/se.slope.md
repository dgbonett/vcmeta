# Computes a slope and its standard error

Computes a slope and its standard error for a simple linear regression
model (random-x model) using the estimated Pearson correlation and the
estimated standard deviations of the response variable and predictor
variable. This function is useful in a meta-analysis of slopes of a
simple linear regression model where some studies report the Pearson
correlation but not the slope.

For more details, see Chapter 1 of Bonett (2021, Volume 5)

## Usage

``` r
se.slope(cor, sdy, sdx, n)
```

## Arguments

- cor:

  estimated Pearson correlation

- sdy:

  estimated standard deviation of the response variable

- sdx:

  estimated standard deviation of the predictor variable

- n:

  sample size

## Value

Returns a one-row matrix:

- Estimate - estimated slope

- SE - standard error

## References

Bonett DG (2021). *Statistical Methods for Psychologists, Vol 1-5,
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
se.slope(.392, 4.54, 2.89, 60)
#>          Estimate        SE
#> Slope:  0.6158062 0.1897647

# Should return: 
#          Estimate        SE
# Slope:  0.6158062 0.1897647

```
