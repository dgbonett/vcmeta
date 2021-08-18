library(vcmeta)

test_that("se.mean2 returns valid matrix", {
  colnames_expected <- c("Estimate", "SE")
  
  res <- se.mean2(21.9, 16.1, 3.82, 3.21, 40, 40)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("se.mean.ps returns valid matrix", {
  colnames_expected <- c("Estimate", "SE")
  
  res <- se.mean.ps(23.9, 25.1, 1.76, 2.01, .78, 25)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("se.stdmean2 returns valid matrix", {
  colnames_expected <- c("Estimate", "SE")
  
  res <- se.stdmean2(21.9, 16.1, 3.82, 3.21, 40, 40, 0)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("se.stdmean.ps returns valid matrix", {
  colnames_expected <- c("Estimate", "SE")
  
  res <- se.stdmean.ps(23.9, 25.1, 1.76, 2.01, .78, 25, 0)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("se.cor returns valid matrix", {
  colnames_expected <- c("Estimate", "SE")
  
  res <- se.cor(.40, 0, 55)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("se.spearman returns valid matrix", {
  colnames_expected <- c("Estimate", "SE")

  res <- se.spear(.40, 55)

  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("se.pbcor returns valid matrix", {
  colnames_expected <- c("Estimate", "SE")
  
  res <- se.pbcor(21.9, 16.1, 3.82, 3.21, 40, 40, 1)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("se.odds returns valid matrix", {
  colnames_expected <- c("Estimate", "SE")
  
  res <- se.odds(36, 50, 21, 50)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("se.meanratio2 returns valid matrix", {
  colnames_expected <- c("Estimate", "SE")
  
  res <- se.meanratio2(21.9, 16.1, 3.82, 3.21, 40, 40)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("se.slope returns valid matrix", {
  colnames_expected <- c("Estimate", "SE")
  
  res <- se.slope(.392, 4.54, 2.89, 60)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})