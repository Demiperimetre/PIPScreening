context("Computation of BF")

test_that("Computation of BF by bridge IS", {
  set.seed(1)
  ech <- runif(100)
  echlogit <- matrix(log(ech/(1+ech)),20,5)
  expect_lt(BFbridgeIS(echlogit,gamma=rep(100,5)), 1)
  expect_equal( computeProbActive(echlogit) ,rep(1,5))
})
