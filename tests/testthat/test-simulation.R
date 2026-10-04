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
