# clonealign

[![Build Status](https://travis-ci.org/kieranrcampbell/clonealign.svg?branch=master)](https://travis-ci.org/kieranrcampbell/clonealign) [![DOI](https://zenodo.org/badge/111455172.svg)](https://zenodo.org/badge/latestdoi/111455172)

`clonealign` assigns single-cell RNA-seq expression to cancer clones by probabilistically mapping RNA-seq to clone-specific copy number profiles using [reparametrization gradient variational inference](https://arxiv.org/abs/1312.6114). This is particularly useful when clones have been inferred using ultra-shallow single-cell DNA-seq meaning SNV analysis is not possible.

<div style="text-align:center">
  <img src="https://raw.githubusercontent.com/kieranrcampbell/clonealign/master/inst/clonealign_figure.png"  align="middle"/>
</div>

## Version 2.0

Clonealign version 2.0 comes with several updated modelling features. In particular:

* A multinomial likelihood that vastly increases runtime and removes the need for custom size factors
* Multiple restarts through the `run_clonealign` function, where the final fit is chosen as that which
maximizes the ELBO

For more info see the NEWS.md file.


## Getting started

### Vignettes

1. [Introduction to clonealign](https://kieranrcampbell.github.io/clonealign/introduction_to_clonealign.html) Overview of `clonealign` including data preparation, model fitting, plotting results, and advanced inference control
2. [Preparing copy number data for input to clonealign](https://kieranrcampbell.github.io/clonealign/preparing_copy_number_data.html) Instructions for taking region/range specific copy number profiles and converting them to gene and clone specific copy numbers for input to clonealign

### Installation

`clonealign` is built using Google's Tensorflow so requires installation of the R package `tensorflow`:

```r
install.packages("tensorflow")
tensorflow::install_tensorflow(extra_packages ="tensorflow-probability", version="1.13.1")
```

Note that `clonealign` uses the [Tensorflow probability](https://www.tensorflow.org/probability/) library, requiring `Tensorflow` version `>= 1.12.0`, which can be installed using the above.

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

[clonealign: statistical integration of independent single-cell RNA and DNA sequencing data from human cancers, _Genome Biology 2019_](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1645-z)

## Authors

Kieran R Campbell, University of British Columbia


