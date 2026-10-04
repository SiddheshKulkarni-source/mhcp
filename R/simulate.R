#' Simulate a Multiplicative Half-Cauchy Process
#'
#' Simulates realizations from the multiplicative half-Cauchy process
#'
#' \deqn{
#' \tilde{\eta}_1 = 1,
#' \qquad
#' \tilde{\eta}_h \sim C^+(0,\zeta), \quad h \ge 2,
#' }
#'
#' with
#'
#' \deqn{
#' \eta_h = \prod_{\ell=1}^h \tilde{\eta}_\ell.
#' }
#'
#' Simulation is performed on the log scale for numerical stability.
#'
#' @param n Number of independent MHCP trajectories.
#' @param d Number of ordered components.
#' @param zeta Positive MHCP scale parameter.
#' @param log Logical. If `TRUE`, return `log(eta)` instead of `eta`.
#'
#' @return An `n` by `d` matrix. Each row is one MHCP trajectory.
#'
#' @examples
#' set.seed(1)
#' x <- rmhcp(n = 5, d = 10, zeta = 0.7)
#' x
#'
#' @export
rmhcp <- function(n, d, zeta, log = FALSE) {

  if (length(n) != 1L || n < 1 || n != as.integer(n)) {
    stop("`n` must be a positive integer.", call. = FALSE)
  }

  if (length(d) != 1L || d < 1 || d != as.integer(d)) {
    stop("`d` must be a positive integer.", call. = FALSE)
  }

  if (length(zeta) != 1L || !is.finite(zeta) || zeta <= 0) {
    stop("`zeta` must be a positive finite number.", call. = FALSE)
  }

  log_eta <- matrix(
    0,
    nrow = n,
    ncol = d
  )

  if (d >= 2) {

    u <- abs(
      matrix(
        stats::rcauchy(n * (d - 1)),
        nrow = n,
        ncol = d - 1
      )
    )

    u <- pmax(
      u,
      .Machine$double.xmin
    )

    increments <- log(zeta) + log(u)

    log_eta[, 2:d] <- t(
      apply(
        increments,
        1,
        cumsum
      )
    )
  }

  colnames(log_eta) <- paste0(
    "eta_",
    seq_len(d)
  )

  if (log) {
    return(log_eta)
  }

  exp(log_eta)
}
