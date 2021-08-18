library(vcmeta)

test_that("meta_ave_mean2 returns valid matrix", {
  colnames_expected <- c("Estimate", "SE", "LL", "UL", "df")
  
  
  m1 <- c(7.4, 6.9)
  m2 <- c(6.3, 5.7)
  sd1 <- c(1.72, 1.53)
  sd2 <- c(2.35, 2.04)
  n1 <- c(40, 60)
  n2 <- c(40, 60)
  res <- meta.ave.mean2(.05, m1, m2, sd1, sd2, n1, n2, bystudy = TRUE)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(3, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("meta.ave.stdmean2 returns valid matrix", {
  colnames_expected <- c("Estimate", "SE", "LL", "UL")
  
  m1 <- c(21.9, 23.1, 19.8)
  m2 <- c(16.1, 17.4, 15.0)
  sd1 <- c(3.82, 3.95, 3.67)
  sd2 <- c(3.21, 3.30, 3.02)
  n1 <- c(40, 30, 24)
  n2 <- c(40, 28, 25)
  res <- meta.ave.stdmean2(.05, m1, m2, sd1, sd2, n1, n2, 0, bystudy = TRUE)
  
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(4, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("meta.ave.mean.ps returns valid matrix", {
  colnames_expected <- c("Estimate", "SE", "LL", "UL", "df")
  
  m1 <- c(53, 60, 53, 57)
  m2 <- c(55, 62, 58, 61)
  sd1 <- c(4.1, 4.2, 4.5, 4.0)
  sd2 <- c(4.2, 4.7, 4.9, 4.8)
  cor <- c(.7, .7, .8, .85)
  n <- c(30, 50, 30, 70)
  res <- meta.ave.mean.ps(.05, m1, m2, sd1, sd2, cor, n, bystudy = TRUE)
  
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(5, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("meta.ave.stdmean.ps returns valid matrix", {
  colnames_expected <- c("Estimate", "SE", "LL", "UL")
  
  m1 <- c(23.9, 24.1)
  m2 <- c(25.1, 26.9)
  sd1 <- c(1.76, 1.58)
  sd2 <- c(2.01, 1.76)
  cor <- c(.78, .84)
  n <- c(25, 30)
  res <- meta.ave.stdmean.ps(.05, m1, m2, sd1, sd2, cor, n, 1, bystudy = TRUE)
  
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(3, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})



test_that("meta.ave.meanratio2 returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "LL", "UL", "exp(Estimate)", "exp(LL)", "exp(UL)", "df"
  )
  
  m1 <- c(7.4, 6.9)
  m2 <- c(6.3, 5.7)
  sd1 <- c(1.7, 1.5)
  sd2 <- c(2.3, 2.0)
  n1 <- c(40, 20)
  n2 <- c(40, 20)
  res <- meta.ave.meanratio2(.05, m1, m2, sd1, sd2, n1, n2, bystudy = TRUE)
  
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(3, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("meta.ave.meanratio.ps returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "LL", "UL", "exp(Estimate)", "exp(LL)", "exp(UL)", "df"
  )
  
  m1 <- c(53, 60, 53, 57)
  m2 <- c(55, 62, 58, 61)
  sd1 <- c(4.1, 4.2, 4.5, 4.0)
  sd2 <- c(4.2, 4.7, 4.9, 4.8)
  cor <- c(.7, .7, .8, .85)
  n <- c(30, 50, 30, 70)
  res <- meta.ave.meanratio.ps(.05, m1, m2, sd1, sd2, cor, n, bystudy = TRUE)
  
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(5, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("meta.ave.cor returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "LL", "UL"
  )
  
  n <- c(55, 190, 65, 35)
  cor <- c(.40, .65, .60, .45)
  res <- meta.ave.cor(.05, n, cor, 0, bystudy = TRUE)
  
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(5, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("meta.ave.slope returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "LL", "UL", "df"
  )
  
  n <- c(45, 85, 50, 60)
  cor <- c(.24, .35, .16, .20)
  sdy <- c(12.2, 14.1, 11.7, 15.9)
  sdx <- c(1.34, 1.87, 2.02, 2.37)
  res <- meta.ave.slope(.05, n, cor, sdy, sdx, bystudy = TRUE)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(5, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("meta.ave.path returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "LL", "UL", "df"
  )
  
  n <- c(75, 85, 250, 160)
  slope <- c(1.57, 1.38, 1.08, 1.25)
  se <- c(.658, .724, .307, .493)
  res <- meta.ave.path(.05, n, slope, se, 2, bystudy = TRUE)

  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(5, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("meta.ave.spear returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "LL", "UL"
  )
  
  n <- c(150, 200, 300, 200, 350)
  cor <- c(.14, .29, .16, .21, .23)
  res <- meta.ave.spear(.05, n, cor, bystudy = TRUE)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(6, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("meta.ave.pbcor returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "LL", "UL"
  )
  
  m1 <- c(21.9, 23.1, 19.8)
  m2 <- c(16.1, 17.4, 15.0)
  sd1 <- c(3.82, 3.95, 3.67)
  sd2 <- c(3.21, 3.30, 3.02)
  n1 <- c(40, 30, 24)
  n2 <- c(40, 28, 25)
  res <- meta.ave.pbcor(.05, m1, m2, sd1, sd2, n1, n2, 2, bystudy = TRUE)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(4, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("meta.ave.semipart returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "LL", "UL"
  )
  
  n <- c(128, 97, 210, 217)
  cor <- c(.35, .41, .44, .39)
  r2 <- c(.29, .33, .36, .39)
  res <- meta.ave.semipart(.05, n, cor, r2, bystudy = TRUE)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(5, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("meta.ave.cronbach returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "LL", "UL"
  )
  
  n <- c(583, 470, 546, 680)
  rel <- c(.91, .89, .90, .89)
  res <- meta.ave.cronbach(.05, n, rel, 10, bystudy = TRUE)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(5, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("meta.ave.odds returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "LL", "UL", "exp(Estimate)", "exp(LL)", "exp(UL)"
  )
  
  n1 <- c(204, 201, 932, 130, 77)
  n2 <- c(106, 103, 415, 132, 83)
  f1 <- c(24, 40, 93, 14, 5)
  f2 <- c(12, 9, 28, 3, 1)
  res <- meta.ave.odds(.05, f1, f2, n1, n2, bystudy = TRUE)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(6, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("meta.ave.propratio2 returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "LL", "UL", "exp(Estimate)", "exp(LL)", "exp(UL)"
  )
  
  n1 <- c(204, 201, 932, 130, 77)
  n2 <- c(106, 103, 415, 132, 83)
  f1 <- c(24, 40, 93, 14, 5)
  f2 <- c(12, 9, 28, 3, 1)
  res <- meta.ave.propratio2(.05, f1, f2, n1, n2, bystudy = TRUE)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(6, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("meta.ave.prop2 returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "LL", "UL"
  )
  
  n1 <- c(204, 201, 932, 130, 77)
  n2 <- c(106, 103, 415, 132, 83)
  f1 <- c(24, 40, 93, 14, 5)
  f2 <- c(12, 9, 28, 3, 1)
  res <- meta.ave.prop2(.05, f1, f2, n1, n2, bystudy = TRUE)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(6, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("meta.ave.prop.ps returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "LL", "UL"
  )
  
  f11 <- c(17, 28, 19)
  f12 <- c(43, 56, 49)
  f21 <- c(3, 5, 5)
  f22 <- c(37, 54, 39)
  res <- meta.ave.prop.ps(.05, f11, f12, f21, f22, bystudy = TRUE)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(4, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("meta.ave.agree returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "LL", "UL"
  )
  
  f11 <- c(17, 28, 19)
  f12 <- c(43, 56, 49)
  f21 <- c(3, 5, 5)
  f22 <- c(37, 54, 39)
  res <- meta.ave.agree(.05, f11, f12, f21, f22, bystudy = TRUE)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(4, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("meta.ave.gen returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "LL", "UL"
  )
  
  est <- c(.022, .751, .421, .287, .052, .146, .562, .904)
  se <- c(.124, .464, .102, .592, .864, .241, .252, .318)
  res <- meta.ave.gen(.05, est, se, bystudy = TRUE)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(9, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("meta.ave.gen.cc returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "LL", "UL"
  )
  
  est <- c(.022, .751, .421, .287, .052, .146, .562, .904)
  se <- c(.124, .464, .102, .592, .864, .241, .252, .318)
  res <- meta.ave.gen.cc(.05, est, se, bystudy = TRUE)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(9, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("meta.ave.gen.rc returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "LL", "UL"
  )
  
  est <- c(.022, .751, .421, .287, .052, .146, .562, .904)
  se <- c(.124, .464, .102, .592, .864, .241, .252, .318)
  res <- meta.ave.gen.rc(.05, est, se, bystudy = TRUE)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(10, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


