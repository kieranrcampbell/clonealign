

#' Plot gene expression and copy number
#'
#' Plot gene expression and copy number as a function of genomic coordinate.
#'
#' This function requires the chromosome, start, and end positions of each features to plot.
#' These are encoded as the \code{chr_str}, \code{start_str}, and \code{end_str} features
#' in \code{rowData(sce)} respectively. If we have ensembl ids or similar we can get add the
#' required fields using the \code{getBMFeatureAnnos} function in the \code{scater} package - \cr
#' \code{
#' sce <- getBMFeatureAnnos(sce, filters = "ensembl_gene_id",
#'                          attributes = c("chromosome_name", "start_position", "end_position"),
#'                          feature_symbol = "hgnc_symbol",
#'                          feature_id = "ensembl_gene_id",
#'                          dataset = "hsapiens_gene_ensembl")
#' } \cr
#' Then we would call \code{plot_clonealign} using \code{chr_str == "chromosome_name"}, \code{start_str == "start_position"},
#' and \code{end_str == "end_position"}.
#'
#' @param sce A \code{SingleCellExperiment}. Note this must have the fields
#' \code{chr_str}, \code{start_str}, and \code{end_str} in \code{rowData}. See details.
#' @param clones The clone assignments of each cell - must be of the same set as the colnames of \code{L}
#' @param cnv_data The gene by clone copy number matrix used as input to \code{clonealign_fit}. The
#' genes (rows) of this must match to the genes in \code{sce}
#' @param chromosome The chromosome to plot
#' @param chr_str The column of \code{rowData(sce)} that refers to the chromosome of each feature
#' @param start_str The column of \code{rowData(sce)} that refers to the start position of each feature
#' @param end_str The column of \code{rowData(sce)} that refers to the end position of each feature
#' @param jitter_cnv Logical - if true random noise is added to the copy number profiles to help
#' display all clones in case two lie on top of each other
#' @param ggplot_palette The palette to be passed to \code{scale_colour_brewer}
#' @param expression_ylim The y axis limits on the smoothed expression plot
#' @param cnv_dodge_sd The standard deviation of the (0 mean) noise to add to the CNV plot
#'
#'
#' @import ggplot2
#' @importFrom stats rnorm sd
#' @importFrom methods is
#' @importFrom dplyr %>%
#' @importFrom SummarizedExperiment rowData "rowData<-"
#'
#' @return A ggplot2 object displaying a CNV and gene expression track for each (inferred) clone
#'
#' @export
#'
#' @examples
#' library(SingleCellExperiment)
#' library(scater)
#' library(dplyr)
#' data(example_clonealign_fit)
#' cnv_data <- rowData(example_sce)[,c("A", "B", "C")]
#' gene_position <- as_data_frame(cnv_data) %>%
#' mutate(gene = seq_len(nrow(cnv_data))) %>%
#'   arrange(A, B, C) %>%
#'  mutate(position = seq_len(nrow(cnv_data))) %>%
#'   arrange(gene) %>%
#'   .$position
#'
#' rowData(example_sce)$chromosome <- "1"
#' rowData(example_sce)$start_pos <- gene_position
#' rowData(example_sce)$end_pos <- gene_position
#' example_sce <- normalize(example_sce)
#' plot_clonealign(example_sce, example_clonealign_fit$clone, cnv_data,
#' chromosome = "1",
#' chr_str = "chromosome",
#' start_str = "start_pos",
#' end_str = "end_pos")
#'
plot_clonealign <- function(sce, clones, cnv_data, chromosome = "1",
                            chr_str = "chr",
                            start_str = "start_position",
                            end_str = "end_position",
                            jitter_cnv = TRUE,
                            ggplot_palette = "Set1",
                            expression_ylim = c(-.15, .15),
                            cnv_dodge_sd = 0.1) {

  scaleFUN <- function(x) sprintf("%.1f", x)

  theme_set(cowplot::theme_cowplot(font_size = 11))

  if(is.data.frame(cnv_data) || is(cnv_data, "DataFrame")) {
    cnv_data <- as.matrix(cnv_data)
  }

  # To appease check() -
  clone <- copy_number <- end <- end_position <- ensembl_gene_id <- NULL
  logcounts <-  mean_exprs <- mean_z_score <- NULL
  per_clone_state_z_score <- position <-  rank_position <- NULL
  sd_exprs <- start <- start_position <- z_score_expression <- NULL

  if(!(chr_str %in% names(rowData(sce)))) {
    stop(paste0("The column 'chr_str' (currently set to '", chr_str, "') must be in rowData(sce) and refer to the chromosome of each gene"))
  }

  if(!(start_str %in% names(rowData(sce)))) {
    stop(paste0("The column 'start_str' (currently set to '", start_str, "') must be in rowData(sce) and refer to the start position of each gene"))
  }

  if(!(end_str %in% names(rowData(sce)))) {
    stop(paste0("The column 'end_str' (currently set to '", end_str, "') must be in rowData(sce) and refer to the end position of each gene"))
  }

  genes_on_chromosome <- rowData(sce)[[ chr_str ]] == chromosome

  if(!any(genes_on_chromosome)) {
    stop(paste("No genes on chromosome", chromosome,  "in CNV regions"))
  }

  sce <- sce[genes_on_chromosome, ]
  cnv_data <- cnv_data[genes_on_chromosome, , drop=FALSE]
  clone_names <- colnames(cnv_data)

  if(!("ensembl_gene_id" %in% names(rowData(sce)))) {
    rowData(sce)$ensembl_gene_id <- as.character(seq_len(nrow(sce)))
  }

  # CNV data
  cnv_df <- as.data.frame(rowData(sce)[,c("ensembl_gene_id", start_str, end_str)])
  cnv_df$rank_position <- rank((cnv_df[[start_str]] + cnv_df[[end_str]])/2)

  cnv_df <- cbind(cnv_df, cnv_data) %>%
    as.data.frame()

  cnv_df[[start_str]] <- cnv_df[[end_str]] <- NULL
  # cnv_df <- tidyr::gather(cnv_df, clone, copy_number, -ensembl_gene_id, -rank_position)


  # cnv_df <- dplyr::mutate(cnv_df, rank_position = rank(position))

  # working out regions

  cnv_spread <- cnv_df# tidyr::spread(cnv_df, clone, copy_number)
  cnv_spread <- dplyr::arrange(cnv_spread, rank_position)
  cnv_spread_clones <- cnv_spread[,clone_names]


  nr <- nrow(cnv_spread_clones)
  state <- 1
  if(nr > 1) {
    for(i in 2:nr) {
      if(all(cnv_spread_clones[i,] == cnv_spread_clones[i-1,])) {
        state <- c(state, state[i-1])
      } else {
        state <- c(state, state[i-1]+1)
      }
    }
  }

  cnv_spread$state <- state

  cnv_df_2_u_grp <- tidyr::gather(cnv_spread, clone, copy_number, -ensembl_gene_id, -rank_position, -state)

  # Need to work out start and end of each region
  cnv_df_2 <- dplyr::group_by(cnv_df_2_u_grp, state, clone, copy_number) %>%
    dplyr::summarise(start = min(rank_position), end = max(rank_position)) %>%
    dplyr::mutate(length = end - start) %>%
    dplyr::ungroup() #%>%
    #dplyr::filter(length > 0)

  if(jitter_cnv) {
    cnv_df_2$copy_number <- cnv_df_2$copy_number + rnorm(length(cnv_df_2$copy_number), 0, cnv_dodge_sd)
  }

  # And that's it for DNA

  dna_plot <- cnv_df_2 %>%
    ggplot(aes(x = start-1, xend = end+1, y = copy_number, yend = copy_number, color = clone)) +
    geom_segment(size = 1.8, position = position_dodge(width = 0.5)) +
    ggplot2::scale_color_brewer(name = "Ground\ntruth\nclone", palette = ggplot_palette) +
    labs(x = "Genomic position", y = "Copy number", subtitle = "scDNA-seq") +
    scale_y_continuous(labels=scaleFUN)

  # RNA data

  lc <- t(as.matrix(logcounts(sce)))
  colnames(lc) <- rowData(sce)$ensembl_gene_id

  gex_df <- dplyr::as_data_frame(lc) %>%
    dplyr::mutate(clone = clones) %>%
    tidyr::gather(ensembl_gene_id, expression, -clone)

  gex_df <- dplyr::inner_join(gex_df, cnv_spread)

  gex_df <- dplyr::filter(gex_df, state %in% cnv_df_2$state)

  mean_sd_df <- dplyr::group_by(gex_df, ensembl_gene_id) %>%
    dplyr::summarise(mean_exprs = mean(expression), sd_exprs = sd(expression))

  gex_df <- dplyr::inner_join(gex_df, mean_sd_df, by = "ensembl_gene_id")

  gex_df$sd_exprs[gex_df$sd_exprs == 0] <- 1

  gex_df <- dplyr::mutate(gex_df, z_score_expression = (expression - mean_exprs) / sd_exprs)

  # gex_df currently normalised z-score values

  norm_gex_per_gene <- dplyr::group_by(gex_df, clone, ensembl_gene_id, rank_position, state) %>%
    dplyr::summarise(mean_z_score = mean(z_score_expression)) %>%
    dplyr::ungroup()

  gex_per_clone_state <- dplyr::group_by(norm_gex_per_gene, clone, state) %>%
    dplyr::summarise(per_clone_state_z_score = mean(mean_z_score)) %>%
    dplyr::ungroup()

  # Now recombine

  df_gex_final <- dplyr::inner_join(norm_gex_per_gene, gex_per_clone_state)


  df_seg <- dplyr::inner_join(cnv_df_2, gex_per_clone_state)

  rna_plot <- df_gex_final %>%
    ggplot(aes(x = rank_position, color = clone)) +
    geom_point(aes(y = mean_z_score), alpha = 0.5) +
    ggplot2::scale_color_brewer(name = "Inferred\nclone", palette = ggplot_palette) +
    labs(x = "Genomic position", y = "Gene expression", subtitle = "scRNA-seq") +
    geom_segment(data = df_seg,
                 aes(x = start-1, xend = end+1, y = per_clone_state_z_score, yend = per_clone_state_z_score),
                 size = 1.2) +
    scale_y_continuous(labels=scaleFUN, limits = expression_ylim)

  cowplot::plot_grid(rna_plot, dna_plot, ncol = 1, align = 'v')

}

