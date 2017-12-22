# clonealign

`clonealign` assigns single-cells to their clones of origin by mapping RNA-seq to clone-specific copy number profiles. This is particularly useful when clones have been inferred using ultra-shallow single-cell DNA-seq meaning SNV analysis is not possible.

<div style="text-align:center">
  <img src="inst/clonealign_figure.png" width="600" align="middle"/>
</div>

## Getting started

### Installation

`clonealign` can be installed from github via

```r
install.packages("devtools") # If not already installed
install_github("kieranrcampbell/clonealign")
```

### Usage

TODO

If `Y` is a cell-by-gene matrix of raw counts and `L` is a gene-by-clone matrix of copy number, then

```r
clonealign(Y, L)
```

## Paper

Coming soon...

## Authors

Kieran R Campbell, University of British Columbia


