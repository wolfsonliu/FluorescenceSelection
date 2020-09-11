# FluorescenceSelection

This R package provides a function `fluorescence.selection` to analyze
fluorescence enriched data.

## Installation

The package can be installed from github by [devtools](https://devtools.r-lib.org/)

```{r}
install.packages("devtools")

devtools::install_github("wolfsonliu/FluorescenceSelection")
```

## Usage

```{r}
result <- fluorescence.selection(
    gene, sgrna, barcode,
    reference, fluorescence,
    selection.ratio=0.05, multiplier=NA, n.to.consider=NA
)
```

`gene`, `sgrna`, `barcode` are character vectors of the library sgRNAs
(with barcode). The `reference` and `fluorescence` are raw counts of
the library, with the same length and order as `gene`, `sgrna`, and
`barcode`.
