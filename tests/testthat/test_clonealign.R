
context("Basic operations")

test_that("clonealign(...) returns a valid object", {
  library(SummarizedExperiment)
  data(example_sce)
  N <- ncol(example_sce)
  G <- nrow(example_sce)
  C <- 3

  L <- rowData(example_sce)[, c("A", "B", "C")]

  cal <- clonealign(example_sce, L, max_iter = 5)

  expect_is(cal, "clonealign_fit")

  clones <- cal$clone

  expect_equal(length(clones), N)

  clonenames_sorted <- sort(unique(clones))

  expect_equal(clonenames_sorted, sort(colnames(L)))

  expect_equal(C, ncol(cal$ml_params$clone_probs))

  expect_equal(N, nrow(cal$ml_params$clone_probs))

  expect_equal(G, length(cal$ml_params$mu))

  expect_true(all(c("mu", "s", "a", "b") %in% names(cal$ml_params)))

})



test_that("clonealign(...) works with non-default size factors", {
  data(example_sce)
  library(SummarizedExperiment)
  N <- ncol(example_sce)

  L <- rowData(example_sce)[, c("A", "B", "C")]

  cal <- clonealign(example_sce, L, max_iter = 5, size_factors = "infer")

  expect_is(cal, "clonealign_fit")

  sf <- colSums(counts(example_sce))

  cal2 <- clonealign(example_sce, L, max_iter = 5, size_factors = sf)
  expect_is(cal2, "clonealign_fit")


})

