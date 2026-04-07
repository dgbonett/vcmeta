# Computes the cell frequencies in a 2x2 table using the marginal proportions and phi correlation

This function computes the cell proportions and frequencies in a 2x2
contingency table using the reported marginal proportions, estimated phi
correlation, and total sample size. The cell frequncies could then be
used to compute other measures of effect size. In the output, "cell ij"
refers to row i and column j.

## Usage

``` r
table.from.phi(p1row, p1col, phi, n)
```

## Arguments

- p1row:

  marginal proportion for row 1

- p1col:

  marginal proportion for column 1

- phi:

  estimated phi correlation

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

## Examples

``` r
table.from.phi(.28, .64, .38, 200)
#>                cell 11   cell 12    cell 21    cell 22
#> Proportion:  0.2610974 0.0189026  0.3789026  0.3410974
#> Frequency:  52.0000000 4.0000000 76.0000000 68.0000000

# Should return:
#                cell 11   cell 12    cell 21    cell 22
# Proportion:  0.2610974 0.0189026  0.3789026  0.3410974
# Frequency   52.0000000 4.0000000 76.0000000 68.0000000

```
