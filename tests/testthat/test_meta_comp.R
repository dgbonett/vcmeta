library(vcmeta)

test_that("meta.sub.cor returns valid matrix", {
  colnames_expected <- c("Estimate", "SE", "LL", "UL")
  
  n <- c(55, 190, 65, 35)
  cor <- c(.40, .65, .60, .45)
  group <- c(1, 1, 2, 0)
  res <- meta.sub.cor(.05, n, cor, 0, group)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(3, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("meta.sub.spear returns valid matrix", {
  colnames_expected <- c("Estimate", "SE", "LL", "UL")
  
  n <- c(55, 190, 65, 35)
  cor <- c(.40, .65, .60, .45)
  group <- c(1, 1, 2, 0)
  res <- meta.sub.spear(.05, n, cor, group)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(3, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("meta.sub.pbcor returns valid matrix", {
  colnames_expected <- c("Estimate", "SE", "LL", "UL")
  
  m1 <- c(45.1, 39.2, 36.3, 34.5)
  m2 <- c(30.0, 35.1, 35.3, 36.2)
  sd1 <- c(10.7, 10.5, 9.4, 11.5)
  sd2 <- c(12.3, 12.0, 10.4, 9.6)
  n1 <- c(40, 20, 50, 25)
  n2 <- c(40, 20, 48, 26)
  group <- c(1, 1, 2, 2)
  res <- meta.sub.pbcor(.05,  m1, m2, sd1, sd2, n1, n2, 2, group)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(3, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("meta.sub.semipart returns valid matrix", {
  colnames_expected <- c("Estimate", "SE", "LL", "UL")
  
  n <- c(55, 190, 65, 35)
  cor <- c(.40, .65, .60, .45)
  r2 <- c(.25, .41, .43, .39)
  group <- c(1, 1, 2, 0)	
  res <- meta.sub.semipart(.05, n, cor, r2, group)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(3, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("meta.sub.cronbach returns valid matrix", {
  colnames_expected <- c("Estimate", "SE", "LL", "UL")
  
  n <- c(120, 170, 150, 135)
  rel <- c(.89, .87, .73, .71)
  group <- c(1, 1, 2, 2)
  q <- 10
  res <- meta.sub.cronbach(.05, n, rel, q, group)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(3, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("meta.lc.mean2 returns valid matrix", {
  colnames_expected <- c("Estimate", "SE", "LL", "UL", "df")
  
  m1 <- c(45.1, 39.2, 36.3, 34.5)
  m2 <- c(30.0, 35.1, 35.3, 36.2)
  sd1 <- c(10.7, 10.5, 9.4, 11.5)
  sd2 <- c(12.3, 12.0, 10.4, 9.6)
  n1 <- c(40, 20, 50, 25)
  n2 <- c(40, 20, 48, 26)
  v <- c(.5, .5, -.5, -.5)
  res <- meta.lc.mean2(.05, m1, m2, sd1, sd2, n1, n2, v)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("meta.lc.stdmean2 returns valid matrix", {
  colnames_expected <- c("Estimate", "SE", "LL", "UL")
  
  m1 <- c(45.1, 39.2, 36.3, 34.5)
  m2 <- c(30.0, 35.1, 35.3, 36.2)
  sd1 <- c(10.7, 10.5, 9.4, 11.5)
  sd2 <- c(12.3, 12.0, 10.4, 9.6)
  n1 <- c(40, 20, 50, 25)
  n2 <- c(40, 20, 48, 26)
  v <- c(.5, .5, -.5, -.5)
  res <- meta.lc.stdmean2(.05, m1, m2, sd1, sd2, n1, n2, v, 0)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("meta.lc.mean.ps returns valid matrix", {
  colnames_expected <- c("Estimate", "SE", "LL", "UL", "df")
  
  m1 <- c(53, 60, 53, 57)
  m2 <- c(55, 62, 58, 61)
  sd1 <- c(4.1, 4.2, 4.5, 4.0)
  sd2 <- c(4.2, 4.7, 4.9, 4.8)
  cor <- c(.7, .7, .8, .85)
  n <- c(30, 50, 30, 70)
  v <- c(.5, .5, -.5, -.5)
  res <- meta.lc.mean.ps(.05, m1, m2, sd1, sd2, cor, n, v)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("meta.lc.meanratio2 returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "LL", "UL", "exp(Estimate)", "exp(LL)", "exp(UL)", "df"
  )
  
  m1 <- c(45.1, 39.2, 36.3, 34.5)
  m2 <- c(30.0, 35.1, 35.3, 36.2)
  sd1 <- c(10.7, 10.5, 9.4, 11.5)
  sd2 <- c(12.3, 12.0, 10.4, 9.6)
  n1 <- c(40, 20, 50, 25)
  n2 <- c(40, 20, 48, 26)
  v <- c(.5, .5, -.5, -.5)
  res <- meta.lc.meanratio2(.05, m1, m2, sd1, sd2, n1, n2, v)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("meta.lc.meanratio.ps returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "LL", "UL", "exp(Estimate)", "exp(LL)", "exp(UL)", "df"
  )
  
  m1 <- c(53, 60, 53, 57)
  m2 <- c(55, 62, 58, 61)
  sd1 <- c(4.1, 4.2, 4.5, 4.0)
  sd2 <- c(4.2, 4.7, 4.9, 4.8)
  cor <- c(.7, .7, .8, .85)
  n <- c(30, 50, 30, 70)
  v <- c(.5, .5, -.5, -.5)
  res <- meta.lc.meanratio.ps(.05, m1, m2, sd1, sd2, cor, n, v)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("meta.lc.odds returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "exp(Estimate)", "exp(LL)", "exp(UL)"
  )
  
  n1 <- c(50, 150, 150)
  f1 <- c(16, 50, 25)
  n2 <- c(50, 150, 150)
  f2 <- c(7, 15, 20)
  v <- c(1, -1, 0)
  res <- meta.lc.odds(.05, f1, f2, n1, n2, v)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("meta.lc.propratio2 returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "exp(Estimate)", "exp(LL)", "exp(UL)"
  )
  
  n1 <- c(50, 150, 150)
  f1 <- c(16, 50, 25)
  n2 <- c(50, 150, 150)
  f2 <- c(7, 15, 20)
  v <- c(1, -1, 0)
  res <- meta.lc.propratio2(.05, f1, f2, n1, n2, v)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("meta.lc.prop2 returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "LL", "UL"
  )
  
  n1 <- c(50, 150, 150)
  f1 <- c(16, 50, 25)
  n2 <- c(50, 150, 150)
  f2 <- c(7, 15, 20)
  v <- c(1, -1, 0)
  res <- meta.lc.prop2(.05, f1, f2, n1, n2, v)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("meta.lc.prop.ps returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "LL", "UL"
  )
  
  f11 <- c(17, 28, 19)
  f12 <- c(43, 56, 49)
  f21 <- c(3, 5, 5)
  f22 <- c(37, 54, 39)
  v <- c(.5, .5, -1)
  res <- meta.lc.prop.ps(.05, f11, f12, f21, f22, v)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("meta.lc.agree returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "LL", "UL"
  )
  
  f11 <- c(17, 28, 19)
  f12 <- c(43, 56, 49)
  f21 <- c(3, 5, 5)
  f22 <- c(37, 54, 39)
  v <- c(.5, .5, -1)
  res <- meta.lc.agree(.05, f11, f12, f21, f22, v)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("meta.lc.mean1 returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "LL", "UL", "df"
  )
  
  m <- c(33.5, 37.9, 38.0, 44.1)
  sd <- c(3.84, 3.84, 3.65, 4.98)
  n <- c(10, 10, 10, 10)
  v <- c(.5, .5, -.5, -.5)
  res <- meta.lc.mean1(.05, m, sd, n, v, eqvar = FALSE)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("meta.lc.prop1 returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "LL", "UL"
  )
  
  f <- c(26, 24, 38)
  n <- c(60, 60, 60)
  v <- c(-.5, -.5, 1)
  res <- meta.lc.prop1(.05, f, n, v)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("meta.lc.gen returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "LL", "UL"
  )
  
  est <- c(.55, .59, .44, .48, .26, .19)
  se <- c(.054, .098, .029, .084, .104, .065)
  v <- c(.5, .5, -.25, -.25, -.25, -.25)
  res <- meta.lc.gen(.05, est, se, v)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})