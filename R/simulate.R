#' Simulate a Multiplicative Half-Cauchy Process
#'
#' Simulates independent trajectories from the multiplicative
#' half-Cauchy process (MHCP).
#'
#' The process is defined by eta_1 = 1 and, for h >= 2,
#' eta_h = eta_(h-1) * tilde_eta_h, where the multiplicative
#' increments follow a half-Cauchy distribution with scale zeta.
#'
#' @param n Number of independent trajectories.
#' @param d Number of ordered components in each trajectory.
#' @param zeta Positive half-Cauchy scale parameter.
#' @param log Logical. If TRUE, return log-scale trajectories.
#'
#' @return A numeric matrix with `n` rows and `d` columns.
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

  if (d >= 2L) {

    u <- abs(
      matrix(
        stats::rcauchy(n * (d - 1L)),
        nrow = n,
        ncol = d - 1L
      )
    )

    u <- pmax(
      u,
      .Machine$double.xmin
    )

    increments <- log(zeta) + log(u)

    log_eta[, 2:d] <- t(
      vapply(
        seq_len(n),
        function(i) {
          cumsum(increments[i, ])
        },
        numeric(d - 1L)
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
