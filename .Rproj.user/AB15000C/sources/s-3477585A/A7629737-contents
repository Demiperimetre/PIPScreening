test_that("output", {
  x <- matrix(runif(80),20,4)
  expect_type(sim3(x,c(.2,.6,.8)), "double")
  expect_equal(length(sim3(x,c(.2,.6,.8))),20)
  expect_type(simGal(x,c(.2,.6,.8)), "double")
  expect_equal(length(simGal(x,c(.2,.6,.8))),20)
})




