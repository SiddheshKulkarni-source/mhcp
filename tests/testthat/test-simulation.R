test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})


test_that("MHCP begins at one", {

  set.seed(1)

  x <- rmhcp(
    n = 100,
    d = 10,
    zeta = 0.7
  )

  expect_true(
    all(
      x[, 1] == 1
    )
  )
})


test_that("simulation returns correct dimensions", {

  x <- rmhcp(
    n = 8,
    d = 12,
    zeta = 0.7
  )

  expect_equal(
    dim(x),
    c(8, 12)
  )
})




test_that("rmhcp handles two-dimensional processes", {

  set.seed(123)

  x <- rmhcp(
    n = 10,
    d = 2,
    zeta = 0.7
  )

  expect_equal(dim(x), c(10, 2))
  expect_true(all(x[, 1] == 1))
  expect_true(all(x[, 2] > 0))
})


test_that("rmhcp handles a single trajectory", {

  set.seed(123)

  x <- rmhcp(
    n = 1,
    d = 10,
    zeta = 0.7
  )

  expect_equal(dim(x), c(1, 10))
  expect_equal(unname(x[1, 1]), 1)
})
