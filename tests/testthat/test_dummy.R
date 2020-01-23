library(HiCDOC)
library(testthat)

context("Dummy check")

test_that("Dummy check", {
    object <- HiCDOCExample()
    expect_is(object, "HiCDOCExp")
})
