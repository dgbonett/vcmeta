library(vcmeta)

test_that("meta.lm.mean2 returns valid matrix", {
  colnames_expected <- c("Estimate", "SE", "t", "p", "LL", "UL", "df")
  
  n1 <- c(65, 30, 29, 45, 50)
  n2 <- c(67, 32, 31, 20, 52)
  m1 <- c(31.1, 32.3, 31.9, 29.7, 33.0)
  m2 <- c(34.1, 33.2, 30.6, 28.7, 26.5)
  sd1 <- c(7.1, 8.1, 7.8, 6.8, 7.6)
  sd2 <- c(7.8, 7.3, 7.5, 7.2, 6.8)
  x1 <- c(4, 6, 7, 7, 8)
  x2 <- c(1, 0, 0, 0, 1)
  X <- matrix(cbind(x1, x2), 5, 2)
  res <- meta.lm.mean2(.05, m1, m2, sd1, sd2, n1, n2, X)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(3, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  testthat::expect_snapshot(res)
})


test_that("meta.lm.stdmean2 returns valid matrix", {
  colnames_expected <- c("Estimate", "SE", "z", "p", "LL", "UL")
  
  n1 <- c(65, 30, 29, 45, 50)
  n2 <- c(67, 32, 31, 20, 52)
  m1 <- c(31.1, 32.3, 31.9, 29.7, 33.0)
  m2 <- c(34.1, 33.2, 30.6, 28.7, 26.5)
  sd1 <- c(7.1, 8.1, 7.8, 6.8, 7.6)
  sd2 <- c(7.8, 7.3, 7.5, 7.2, 6.8)
  x1 <- c(4, 6, 7, 7, 8)
  X <- matrix(x1, 5, 1)
  res <- meta.lm.stdmean2(.05, m1, m2, sd1, sd2, n1, n2, X, 0)
  
  
  n1 <- c(65, 30, 29, 45, 50)
  n2 <- c(67, 32, 31, 20, 52)
  m1 <- c(31.1, 32.3, 31.9, 29.7, 33.0)
  m2 <- c(34.1, 33.2, 30.6, 28.7, 26.5)
  sd1 <- c(7.1, 8.1, 7.8, 6.8, 7.6)
  sd2 <- c(7.8, 7.3, 7.5, 7.2, 6.8)
  x1 <- c(4, 6, 7, 7, 8)
  X <- matrix(x1, 5, 1)
  meta.lm.stdmean2(.05, m1, m2, sd1, sd2, n1, n2, X, 0)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(2, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  testthat::expect_snapshot(res)
})



test_that("meta.lm.mean.ps returns valid matrix", {
  colnames_expected <- c("Estimate", "SE", "t", "p", "LL", "UL", "df")
  
  n <- c(65, 30, 29, 45, 50)
  cor <- c(.87, .92, .85, .90, .88)
  m1 <- c(20.1, 20.5, 19.3, 21.5, 19.4)
  m2 <- c(10.4, 10.2, 8.5, 10.3, 7.8)
  sd1 <- c(9.3, 9.9, 10.1, 10.5, 9.8)
  sd2 <- c(7.8, 8.0, 8.4, 8.1, 8.7)
  x1 <- c(2, 3, 3, 4, 4)
  X <- matrix(x1, 5, 1)
  res <- meta.lm.mean.ps(.05, m1, m2, sd1, sd2, cor, n, X)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(2, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  testthat::expect_snapshot(res)
})


test_that("meta.lm.stdmean.ps returns valid matrix", {
  colnames_expected <- c("Estimate", "SE", "z", "p", "LL", "UL")
  
  n <- c(65, 30, 29, 45, 50)
  cor <- c(.87, .92, .85, .90, .88)
  m1 <- c(20.1, 20.5, 19.3, 21.5, 19.4)
  m2 <- c(10.4, 10.2, 8.5, 10.3, 7.8)
  sd1 <- c(9.3, 9.9, 10.1, 10.5, 9.8)
  sd2 <- c(7.8, 8.0, 8.4, 8.1, 8.7)
  x1 <- c(2, 3, 3, 4, 4)
  X <- matrix(x1, 5, 1)
  res <- meta.lm.stdmean.ps(.05, m1, m2, sd1, sd2, cor, n, X, 0)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(2, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  testthat::expect_snapshot(res)
})


test_that("meta.lm.meanratio2 returns valid matrix", {
  colnames_expected <- c(
    "Estimate", 
    "SE", 
    "z", 
    "p", 
    "LL", 
    "UL", 
    "exp(Estimate)", 
    "exp(LL)", 
    "exp(UL)"
  )
  
  n1 <- c(65, 30, 29, 45, 50)
  n2 <- c(67, 32, 31, 20, 52)
  m1 <- c(31.1, 32.3, 31.9, 29.7, 33.0)
  m2 <- c(34.1, 33.2, 30.6, 28.7, 26.5)
  sd1 <- c(7.1, 8.1, 7.8, 6.8, 7.6)
  sd2 <- c(7.8, 7.3, 7.5, 7.2, 6.8)
  x1 <- c(4, 6, 7, 7, 8)
  X <- matrix(x1, 5, 1)
  res <- meta.lm.meanratio2(.05, m1, m2, sd1, sd2, n1, n2, X)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(2, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  testthat::expect_snapshot(res)
})


test_that("meta.lm.meanratio.ps returns valid matrix", {
  colnames_expected <- c(
    "Estimate", 
    "SE", 
    "z", 
    "p", 
    "LL", 
    "UL", 
    "exp(Estimate)", 
    "exp(LL)", 
    "exp(UL)"
  )
  
  n <- c(65, 30, 29, 45, 50)
  cor <- c(.87, .92, .85, .90, .88)
  m1 <- c(20.1, 20.5, 19.3, 21.5, 19.4)
  m2 <- c(10.4, 10.2, 8.5, 10.3, 7.8)
  sd1 <- c(9.3, 9.9, 10.1, 10.5, 9.8)
  sd2 <- c(7.8, 8.0, 8.4, 8.1, 8.7)
  x1 <- c(2, 3, 3, 4, 4)
  X <- matrix(x1, 5, 1)
  res <- meta.lm.meanratio.ps(.05, m1, m2, sd1, sd2, cor, n, X)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(2, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  testthat::expect_snapshot(res)
})


test_that("meta.lm.cor.gen returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "z", "p"
  )
  
  cor <- c(.40, .65, .60, .45)
  se <- c(.182, .114, .098, .132)
  x1 <- c(18, 25, 23, 19)
  X <- matrix(x1, 4, 1)
  res <- meta.lm.cor.gen(.05, cor, se, X)
  
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(2, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  testthat::expect_snapshot(res)
})


test_that("meta.lm.cor returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "z", "p", "LL", "UL"
  )
  
  n <- c(55, 190, 65, 35)
  cor <- c(.40, .65, .60, .45)
  q <- 0
  x1 <- c(18, 25, 23, 19)
  X <- matrix(x1, 4, 1)
  res <- meta.lm.cor(.05, n, cor, q, X)
  
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(2, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  testthat::expect_snapshot(res)
})


test_that("meta.lm.spear returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "z", "p", "LL", "UL"
  )
  
  n <- c(150, 200, 300, 200, 350)
  cor <- c(.14, .29, .16, .21, .23)
  x1 <- c(18, 25, 23, 19, 24)
  X <- matrix(x1, 5, 1)
  res <- meta.lm.spear(.05, n, cor, X)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(2, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  testthat::expect_snapshot(res)
})


test_that("meta.lm.spear returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "z", "p", "LL", "UL"
  )
  
  n <- c(150, 200, 300, 200, 350)
  cor <- c(.14, .29, .16, .21, .23)
  x1 <- c(18, 25, 23, 19, 24)
  X <- matrix(x1, 5, 1)
  res <- meta.lm.spear(.05, n, cor, X)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(2, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  testthat::expect_snapshot(res)
})


test_that("meta.lm.semipart returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "z", "p", "LL", "UL"
  )
  
  n <- c(128, 97, 210, 217)
  cor <- c(.35, .41, .44, .39)
  r2 <- c(.29, .33, .36, .39)
  x1 <- c(18, 25, 23, 19)
  X <- matrix(x1, 4, 1)
  res <- meta.lm.semipart(.05, n, cor, r2, X)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(2, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  testthat::expect_snapshot(res)
})


test_that("meta.lm.cronbach returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "z", "p", "LL", "UL"
  )
  
  n <- c(583, 470, 546, 680)
  rel <- c(.91, .89, .90, .89)
  x1 <- c(1, 0, 0, 0)
  X <- matrix(x1, 4, 1)
  res <- meta.lm.cronbach(.05, n, rel, 10, X)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(2, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  testthat::expect_snapshot(res)
})


test_that("meta.lm.odds returns valid matrix", {
  colnames_expected <- c(
    "Estimate", 
    "SE", 
    "z", 
    "p", 
    "LL", 
    "UL", 
    "exp(Estimate)", 
    "exp(LL)", 
    "exp(UL)"
  )
  
  n1 <- c(204, 201, 932, 130, 77)
  n2 <- c(106, 103, 415, 132, 83)
  f1 <- c(24, 40, 93, 14, 5)
  f2 <- c(12, 9, 28, 3, 1)
  x1 <- c(4, 4, 5, 3, 26)
  x2 <- c(1, 1, 1, 0, 0)
  X <- matrix(cbind(x1, x2), 5, 2)
  res <- meta.lm.oddsratio(.05, f1, f2, n1, n2, X)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(3, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("meta.lm.propratio2 returns valid matrix", {
  colnames_expected <- c(
    "Estimate", 
    "SE", 
    "z", 
    "p", 
    "LL", 
    "UL", 
    "exp(Estimate)", 
    "exp(LL)", 
    "exp(UL)"
  )
  
  n1 <- c(204, 201, 932, 130, 77)
  n2 <- c(106, 103, 415, 132, 83)
  f1 <- c(24, 40, 93, 14, 5)
  f2 <- c(12, 9, 28, 3, 1)
  x1 <- c(4, 4, 5, 3, 26)
  x2 <- c(1, 1, 1, 0, 0)
  X <- matrix(cbind(x1, x2), 5, 2)
  res <- meta.lm.propratio2(.05, f1, f2, n1, n2, X)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(3, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  testthat::expect_snapshot(res)
})


test_that("meta.lm.prop2 returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "z", "p", "LL", "UL"
  )
  
  f1 <- c(24, 40, 93, 14, 5)
  f2 <- c(12, 9, 28, 3, 1)
  n1 <- c(204, 201, 932, 130, 77)
  n2 <- c(106, 103, 415, 132, 83)
  x1 <- c(4, 4, 5, 3, 26)
  x2 <- c(1, 1, 1, 0, 0)
  X <- matrix(cbind(x1, x2), 5, 2)
  res <- meta.lm.prop2(.05, f1, f2, n1, n2, X)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(3, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  testthat::expect_snapshot(res)
})


test_that("meta.lm.prop.ps returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "z", "p", "LL", "UL"
  )
  
  f11 <- c(40, 20, 25, 30)
  f12 <- c(3, 2, 2, 1)
  f21 <- c(7, 6, 8, 6)
  f22 <- c(26, 25, 13, 25)
  x1 <- c(1, 1, 4, 6)
  x2 <- c(1, 1, 0, 0)
  X <- matrix(cbind(x1, x2), 4, 2)
  res <- meta.lm.prop.ps(.05, f11, f12, f21, f22, X)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(3, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  testthat::expect_snapshot(res)
})


test_that("meta.lm.agree returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "z", "p", "LL", "UL"
  )
  
  f11 <- c(40, 20, 25, 30)
  f12 <- c(3, 2, 2, 1)
  f21 <- c(7, 6, 8, 6)
  f22 <- c(26, 25, 13, 25)
  x1 <- c(1, 1, 4, 6)
  x2 <- c(1, 1, 0, 0)
  X <- matrix(cbind(x1, x2), 4, 2)
  res <- meta.lm.agree(.05, f11, f12, f21, f22, X)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(3, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  testthat::expect_snapshot(res)
})


test_that("meta.lm.mean1 returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "t", "p", "LL", "UL", "df"
  )
  
  n <- c(25, 15, 30, 25, 40)
  m <- c(20.1, 20.5, 19.3, 21.5, 19.4)
  sd <- c(10.4, 10.2, 8.5, 10.3, 7.8)
  x1 <- c(1, 1, 0, 0, 0)
  x2 <- c( 12, 13, 11, 13, 15)
  X <- matrix(cbind(x1, x2), 5, 2)
  res <- meta.lm.mean1(.05, m, sd, n, X)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(3, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  testthat::expect_snapshot(res)
})


test_that("meta.lm.prop1 returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "z", "p", "LL", "UL"
  )
  
  f <- c(38, 26, 24, 15, 45, 38)
  n <- c(80, 60, 70, 50, 180, 200)
  x1 <- c(10, 15, 18, 22, 24, 30)
  X <- matrix(x1, 6, 1)
  res <- meta.lm.prop1(.05, f, n, X)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(2, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  testthat::expect_snapshot(res)
})


test_that("meta.lm.gen returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "z", "p", "LL", "UL"
  )
  
  est <- c(4.1, 4.7, 4.9, 5.7, 6.6, 7.3)
  se <- c(1.2, 1.5, 1.3, 1.8, 2.0, 2.6)
  x1 <- c(10, 20, 30, 40, 50, 60)
  x2 <- c(1, 1, 1, 0, 0, 0)
  X <- matrix(cbind(x1, x2), 6, 2)
  res <- meta.lm.gen(.05, est, se, X)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(3, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  testthat::expect_snapshot(res)
})
