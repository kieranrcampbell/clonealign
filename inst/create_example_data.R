library(scater)

sce <- readRDS("../data/scrnaseq/SA501X2B/SA501X2B_10x_cnv_no_X_use_nov_30.rds")

set.seed(123)
cells_to_sample <- sample(seq_len(ncol(sce)), 200)
genes_to_sample <- sample(seq_len(nrow(sce)), 100)
sce_example <- sce[genes_to_sample, cells_to_sample]

Y <- as.matrix(counts(sce_example))

rownames(Y) <- paste0("gene_", seq_len(nrow(Y)))
colnames(Y) <- paste0("cell_", seq_len(ncol(Y)))

L <- rowData(sce_example)[, c("A", "B", "C")]
rownames(L) <- rownames(Y)

example_sce <- SingleCellExperiment(assays = list(counts = Y),
                          rowData = L)

use_data(example_sce, pkg = ".")
