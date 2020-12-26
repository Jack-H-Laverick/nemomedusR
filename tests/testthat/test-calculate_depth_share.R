test_that("Shares are correct using default boundaries", {
  expect_equal(calculate_depth_share(seq(0, 100, by = 10)), 
               c(5, rep(10, 9), 5))
})
