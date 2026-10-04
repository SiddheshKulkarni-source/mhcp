#' Expected Threshold-Based Effective Dimension
#'
#' Computes the expected number of MHCP components exceeding
#' a threshold `epsilon`,
#'
#' \deqn{
#' E(N_\epsilon)
#' =
#' \sum_{h=1}^{\infty}
#' P(\eta_h > \epsilon).
#' }
#'
#' The infinite sum is numerically truncated once exceedance
#' probabilities become negligible.
#'
#' @param zeta Positive MHCP parameter.
#' @param epsilon Positive threshold.
#' @param max_h Maximum component index considered.
#' @param tol Numerical stopping tolerance.
#' @param consecutive Number of consecutive terms below `tol`
#'   required before stopping.
#'
#' @return A list containing the estimated expected dimension,
#'   the truncation index, and the individual exceedance probabilities.
#'
#' @examples
#' expected_mhcp_dimension(
#'   zeta = 0.369411,
#'   epsilon = 0.05
#' )
#'
#' @export
expected_mhcp_dimension <- function(
    zeta,
    epsilon,
    max_h = 500,
    tol = 1e-10,
    consecutive = 10
) {

  if (zeta <= 0 || zeta >= 1) {
    stop(
      "`zeta` must satisfy 0 < zeta < 1 for the shrinking regime.",
      call. = FALSE
    )
  }

  if (epsilon <= 0) {
    stop("`epsilon` must be positive.", call. = FALSE)
  }

  if (max_h < 2) {
    stop("`max_h` must be at least 2.", call. = FALSE)
  }

  exceedance <- numeric(max_h)

  # eta_1 = 1 deterministically
  exceedance[1] <- as.numeric(
    1 > epsilon
  )

  small_count <- 0L
  stopping_h <- max_h

  for (h in 2:max_h) {

    prob <- 1 -
      pmhcp(
        x = epsilon,
        H = h,
        zeta = zeta
      )

    exceedance[h] <- prob

    if (prob < tol) {
      small_count <- small_count + 1L
    } else {
      small_count <- 0L
    }

    if (small_count >= consecutive) {
      stopping_h <- h
      break
    }
  }

  exceedance <- exceedance[
    seq_len(stopping_h)
  ]

  list(
    mean_effective_dimension =
      sum(exceedance),

    last_index =
      stopping_h,

    threshold =
      epsilon,

    zeta =
      zeta,

    exceedance_probability =
      exceedance
  )
}
