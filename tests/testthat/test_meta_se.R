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

test_that("se.prop2 returns valid matrix", {
  colnames_expected <- c("Estimate", "SE")
  
  res <- se.prop2(31, 16, 40, 40)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})

test_that("se.prop.ps returns valid matrix", {
  colnames_expected <- c("Estimate", "SE")

  res <- se.prop.ps(16, 64, 5, 15)

  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("se.ave.mean2.dep returns valid matrix", {
  colnames_expected <- c("Estimate",        "SE",    "VAR(A)",    "VAR(B)",  "COV(A,B)")
  
  res <- se.ave.mean2.dep(21.9, 16.1, 3.82, 3.21, 24.8, 17.1, 3.57, 3.64, .785, 40, 40)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("se.ave.cor.over returns valid matrix", {
  colnames_expected <- c("Estimate",        "SE",   "VAR(cor12)", "VAR(cor13)", "COV(cor12,cor13)")
  
  res <- se.ave.cor.over(.462, .518, .755, 100)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(2, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("se.ave.cor.nonover returns valid matrix", {
  colnames_expected <- c("Estimate",        "SE",   "VAR(cor12)", "VAR(cor34)", "COV(cor12,cor34)")
  
  res <- se.ave.cor.nonover(.357, .398, .755, .331, .347, .821, 100)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(2, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("se.tetra returns valid matrix", {
  colnames_expected <- c("Estimate",        "SE")
  
  res <- se.tetra(46, 15, 54, 85)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("se.biphi returns valid matrix", {
  colnames_expected <- c("Estimate",        "SE")
  
  res <- se.biphi(34, 22, 50, 50)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("se.cohen returns valid matrix", {
  colnames_expected <- c("Estimate",        "SE")
  
  res <- se.cohen(.78, 35, 50)

  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})