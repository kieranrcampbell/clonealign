
#' Check and pre-process gene expression data
#' @keywords internal
#' 
#' @return The processed gene expression data
process_gene_expression_data <- function(gene_expression_data) {
  Y <- NULL
  
  # Parse gene expression data first
  if(is(gene_expression_data, "SingleCellExperiment") || is(gene_expression_data, "SummarizedExperiment")) {
    assay_names <- names(assays(gene_expression_data))
    if(!("counts" %in% assay_names)) {
      stop(paste("counts not in assays(gene_expression_data). Available assays:", paste(assay_names, collapse = ",")))
    }
    Y <- t(as.matrix(assay(gene_expression_data, "counts")))
  } else if(is(gene_expression_data, "dgCMatrix")) {
    Y <- as.matrix(gene_expression_data)
  } else if(is(gene_expression_data, "matrix")) {
    Y <- gene_expression_data
  } else {
    stop("Input gene_expression_data must be SingleCellExperiment, SummarizedExperiment, or matrix")
  }
  
  Y
}

#' Preprocess copy number data
#' @keywords internal
#' 
#' @return The processed copy number data
process_cnv_data <- function(copy_number_data) {
  L <- NULL
  if(is(copy_number_data, "data.frame") || is(copy_number_data, "DataFrame")) {
    L <- as.matrix(copy_number_data)
  } else if(is(copy_number_data, "matrix")) {
    L <- copy_number_data
  } else {
    stop(paste("copy_number_data must be a matrix, data.frame or DataFrame. Current class:", class(copy_number_data)))
  }
  L
}

#' Get outlier genes
#' 
#' @keywords internal
#' 
#' @param Y Cell by gene matrix of counts
#' @param nmads Number of mads above which gene is considered outlier
#' 
#' @return A logical vector whether each gene (as represented by a column of Y) 
#' is an outlier 
#' @importFrom stats mad
#' @importFrom matrixStats rowMaxs
#' 
#' @examples 
#' Y <- matrix(rpois(200, 10), 20, 10)
#' clonealign:::get_outlying_genes(Y, 3)
get_outlying_genes <- function(Y, nmads) {
  gene_means <- colMeans(Y)
  md <- mad(gene_means)
  gene_means > mean(gene_means) + nmads * md
}

#' Preprocess data for input to clonealign
#' 
#' Preprocess data for input to clonealign, filtering for cells and genes with minimum counts, along
#' with removing outlying genes, and genes with the same copy number between clones.
#' 
#' @param gene_expression_data Input gene expression data. See \code{clonealign} for details
#' @param copy_number_data Input copy number data. See \code{clonealign} for details
#' @param min_counts_per_gene Minimum counts per gene for the gene to pass filtering
#' @param min_counts_per_cell Minimum counts per cell for the cell to pass filtering
#' @param remove_outlying_genes Logical - should genes whose expression is an outlier wrt
#' all others be removed?
#' @param nmads The number of median absolute deviations (MADs) the per-gene mean of the 
#' raw counts is from the overall mean to be considered an outlier 
#' @param remove_genes_same_copy_number Logical - should genes with the same copy number in all clones be removed?
#' @param max_copy_number Maximum copy number per gene to retain (see "saturation" under original paper)
#' 
#' @export
#' 
#' @return A list with entries \code{gene_expression_data} and \code{copy_number_data}
#' for input to clonealign, along with names of the retained genes and cells after filtering.
#' 
#' @importFrom matrixStats rowMaxs
#' @importFrom matrixStats rowVars
#' 
#' @examples 
#' library(SummarizedExperiment)
#' data("example_sce")
#' L <- rowData(example_sce)[,c("A", "B", "C")]
#' ca_data <- preprocess_for_clonealign(example_sce, L)
preprocess_for_clonealign <- function(gene_expression_data,
                                  copy_number_data,
                                  min_counts_per_gene = 20,
                                  min_counts_per_cell = 100,
                                  remove_outlying_genes = TRUE,
                                  nmads = 10,
                                  max_copy_number = 6,
                                  remove_genes_same_copy_number = TRUE) {
  
 
  Y <- process_gene_expression_data(gene_expression_data)
  G <- ncol(Y)
  
  L <- process_cnv_data(copy_number_data)
  
  if(nrow(L) != G) {
    stop("copy_number_data must have same number of genes (rows) as gene_expression_data")
  }
  

  
  copy_number_exceeds_max <- rowMaxs(L) > max_copy_number
  Y <- Y[, !copy_number_exceeds_max]
  L <- L[!copy_number_exceeds_max,]
  
  expressed_sufficiently <- colSums(Y) > min_counts_per_gene
  Y <- Y[, expressed_sufficiently]
  L <- L[expressed_sufficiently, ]
  
  # Remove outlier genes
  if(remove_outlying_genes) {
    outlying_genes <- get_outlying_genes(Y, nmads)

    Y <- Y[, !outlying_genes]
    L <- L[!outlying_genes,]
  }
  
  ## Remove genes with same copy number
  if(remove_genes_same_copy_number) {
    to_remove <- rowVars(L) == 0
    Y <- Y[, !to_remove]
    L <- L[!to_remove, ]
  }
    
  ## Remove cells with no counts mapping
  cells_with_coverage <- rowSums(Y) > min_counts_per_cell
  Y <- Y[cells_with_coverage, ]
  
  list(
    gene_expression_data = Y,
    copy_number_data = L,
    retained_cells = rownames(Y),
    retained_genes = colnames(Y)
  )
}