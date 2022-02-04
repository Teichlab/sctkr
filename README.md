
<!-- README.md is generated from README.Rmd. Please edit that file -->
sctkr
=====

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental) <!-- [![CRAN status](https://www.r-pkg.org/badges/version/sctkr)](https://CRAN.R-project.org/package=sctkr) --> <!-- badges: end -->

The goal of sctkr is to serve as a deposit for functions that perform common tasks for single cell analysis

Installation
------------

You can install the released version of sctkr from [github](https://github.com/Teichlab/sctkr) with:

``` r
devtools::install_github("Teichlab/sctkr")
```

Usage example
------------

A simple usage example is given below, for more options, look at "?CellTypeCompositionAnalysis" and "?plot_ranef"

``` r
result <- CellTypeCompositionAnalysis(
    obs_tbl,
    'sample_id',
    'my_annot',
    colVarCats=c(
        'Sex', 'Kit_version', 'Sample_location', 'Age_bin'
    )
)

print(plot_ranef(
    result$ranef,
    vars=list(
        Sample_location=c('Nose', 'Trachea', 'Bronchi'),
        Age_bin=c('Neonate', 'Infant', 'Young child', 'Child', 'Adolescent', 'Adult', 'Eldly')
    ),
    maxFC=3
))
```
