test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})


test_that("centered process has median zero", {

  expect_equal(
    pmhcp_centered(
      0,
      H = 10
    ),
    0.5,
    tolerance = 1e-10
  )
})


test_that("90 percent calibration reproduces target probability", {

  result <- calibrate_mhcp(
    H = 10,
    epsilon = 0.05,
    p = 0.90
  )

  probability <- pmhcp(
    x = 0.05,
    H = 10,
    zeta = result$zeta
  )

  expect_equal(
    probability,
    0.90,
    tolerance = 1e-6
  )
})


test_that("median calibration has closed form", {

  result <- calibrate_mhcp(
    H = 10,
    epsilon = 0.05,
    p = 0.5
  )

  theoretical <-
    0.05^(1 / 9)

  expect_equal(
    result$zeta,
    theoretical,
    tolerance = 1e-8
  )
})


test_that("worked example reproduces manuscript value", {

  result <- calibrate_mhcp(
    H = 10,
    epsilon = 0.05,
    p = 0.90
  )

  expect_equal(
    result$zeta,
    0.3694,
    tolerance = 5e-4
  )
})
