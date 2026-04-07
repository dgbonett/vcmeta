# Computes the cell frequencies in a 2x2 table using the marginal proportions and odds ratio

This function computes the cell proportions and frequencies in a 2x2
contingency table using the reported marginal proportions, estimated
odds ratio, and total sample size. The cell frequncies could then be
used to compute other measures of effect size. In the output, "cell ij"
refers to row i and column j.

## Usage

``` r
table.from.odds(p1row, p1col, or, n)
```

## Arguments

- p1row:

  marginal proportion for row 1

- p1col:

  marginal proportion for column 1

- or:

  estimated odds ratio

- n:

  total sample size

## Value

A 2-row matrix. The rows are:

- Row 1 gives the four computed cell proportions

- Row 2 gives the four computed cell frequencies

The columns are:

- cell 11 - proportion and frequency for cell 11

- cell 12 - proportion and frequency for cell 12

- cell 21 - proportion and frequency for cell 21

- cell 22 - proportion and frequency for cell 22

## References

Bonett DG (2007). “Transforming odds ratios into correlations for
meta-analytic research.” *American Psychologist*, **62**(3), 254–255.
[doi:10.1037/0003-066X.62.3.254](https://doi.org/10.1037/0003-066X.62.3.254)
.

## Examples

``` r
table.from.odds(.17, .5, 3.18, 100)
#>                cell 11    cell 12    cell 21    cell 22
#> Proportion:  0.1233262 0.04667383  0.3766738  0.4533262
#> Frequency:  12.0000000 5.00000000 38.0000000 45.0000000

# Should return:
#                cell 11    cell 12    cell 21    cell 22
# Proportion:  0.1233262 0.04667383  0.3766738  0.4533262
# Frequency:  12.0000000 5.00000000 38.0000000 45.0000000

```
