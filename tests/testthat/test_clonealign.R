
context("Basic operations")

test_that("clonealign(...) returns a valid object", {
  library(SummarizedExperiment)
  data(example_sce)
  N <- ncol(example_sce)
  G <- nrow(example_sce)
  C <- 3

  L <- rowData(example_sce)[, c("A", "B", "C")]

  suppressWarnings({
    cal <- clonealign(example_sce, L, max_iter = 5)
  })

  expect_is(cal, "clonealign_fit")

  clones <- cal$clone

  expect_equal(length(clones), N)

  inferred_clones <- unique(clones)

  expect_true(all(inferred_clones %in% c(colnames(L), "unassigned")))

  expect_equal(C, ncol(cal$ml_params$clone_probs))

  expect_equal(N, nrow(cal$ml_params$clone_probs))

  expect_equal(length(cal$retained_genes), length(cal$ml_params$mu))
  
  expect_lte(length(cal$ml_params$mu), G)

  expect_true(all(c("clone_probs", "mu", "s") %in% names(cal$ml_params)))
  
  expect_true(all(c("clone", "convergence_info", "retained_genes", "correlations", "ml_params") %in% names(cal)))

})


test_that("Seed setting works correctly", {
  library(SummarizedExperiment)
  data(example_sce)
  N <- ncol(example_sce)
  G <- nrow(example_sce)
  C <- 3
  
  L <- rowData(example_sce)[, c("A", "B", "C")]
  
  suppressWarnings({
    set.seed(12345L)
    cal1 <- clonealign(example_sce, L, max_iter = 5)
  })
  
  suppressWarnings({
    set.seed(12345L)
    cal2 <- clonealign(example_sce, L, max_iter = 5)
  })
  
  expect_equal(
    cal1$convergence_info$final_elbo,
    cal2$convergence_info$final_elbo
  )

})
