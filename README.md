# clonealign

[![Build Status](https://travis-ci.org/kieranrcampbell/clonealign.svg?branch=master)](https://travis-ci.org/kieranrcampbell/clonealign)

`clonealign` assigns single-cell RNA-seq expression to cancer clones by mapping RNA-seq to clone-specific copy number profiles. This is particularly useful when clones have been inferred using ultra-shallow single-cell DNA-seq meaning SNV analysis is not possible.

<div style="text-align:center">
  <img src="inst/clonealign_idea.png" width="500" align="middle"/>
</div>

## Getting started

### Installation

`clonealign` is built using Google's Tensorflow so requires installation of the R package `tensorflow`:

```r
install.packages("tensorflow")
tensorflow::install_tensorflow()
```

`clonealign` can then be installed from github:

```r
install.packages("devtools") # If not already installed
install_github("kieranrcampbell/clonealign")
```

### Usage

`clonealign` accepts either a cell-by-gene matrix of raw counts or a [SingleCellExperiment](https://bioconductor.org/packages/3.7/bioc/html/SingleCellExperiment.html) with a `counts` assay as gene expression input. It also requires a gene-by-clone matrix or `data.frame` corresponding to the copy number of each gene in each clone. The cells are then assigned to their clones by calling

```r
cal <- clonealign(gene_expression_data, # matrix or SingleCellExperiment
                  copy_number_data)     # matrix or data.frame
print(cal)
```
```
A clonealign_fit for 200 cells, 100 genes, and 3 clones
To access clone assignments, call x$clone
To access ML parameter estimates, call x$ml_params
```

```r
print(head(cal$clone))
```
```
[1] "B" "C" "C" "B" "C" "B"
```


## Paper

Coming soon...


## Authors

Kieran R Campbell, University of British Columbia


