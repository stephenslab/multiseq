test_that("default ash parameters are correct", {
  expect_equal(setAshParam(list()),readRDS("default.ashparam.RDS"))
})