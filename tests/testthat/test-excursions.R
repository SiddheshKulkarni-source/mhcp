test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})


test_that("excursion probabilities agree with theory", {

  p <- mhcp_excursion_prob(
    zeta = 0.7,
    c = 1
  )

  expect_equal(
    p,
    (2 / pi) * atan(0.7),
    tolerance = 1e-12
  )
})
