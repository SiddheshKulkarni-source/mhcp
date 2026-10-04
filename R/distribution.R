#' CDF of the Centered MHCP Log Process
#'
#' Evaluates the distribution function of
#'
#' \deqn{
#' Z_H = \sum_{\ell=2}^{H} Y_\ell,
#' }
#'
#' where
#'
#' \deqn{
#' Y_\ell = \log U_\ell,
#' \qquad
#' U_\ell \sim C^+(0,1).
#' }
#'
#' The characteristic function is
#'
#' \deqn{
#' \phi_{Z_H}(t)
#' =
#' \operatorname{sech}^{H-1}
#' \left(\frac{\pi t}{2}\right).
#' }
#'
#' The CDF is evaluated using Gil-Pelaez inversion.
#'
#' @param q Numeric vector of evaluation points.
#' @param H Component index. Must be an integer at least 2.
#' @param rel_tol Relative tolerance used by numerical integration.
#' @param upper Upper integration limit.
#'
#' @return Numeric vector containing `P(Z_H <= q)`.
#'
#' @examples
#' pmhcp_centered(0, H = 10)
#'
#' @export
pmhcp_centered <- function(
    q,
    H,
    rel_tol = 1e-8,
    upper = 30
) {

  if (length(H) != 1L || H < 2 || H != as.integer(H)) {
    stop("`H` must be an integer >= 2.", call. = FALSE)
  }

  if (rel_tol <= 0) {
    stop("`rel_tol` must be positive.", call. = FALSE)
  }

  n_inc <- H - 1

  one_cdf <- function(z) {

    if (is.na(z)) {
      return(NA_real_)
    }

    if (is.infinite(z)) {
      return(if (z < 0) 0 else 1)
    }

    if (abs(z) < 1e-14) {
      return(0.5)
    }

    integrand <- function(t) {

      sinc_term <- ifelse(
        abs(t) < 1e-12,
        z,
        sin(t * z) / t
      )

      sinc_term *
        sech(pi * t / 2)^n_inc /
        pi
    }

    integral_value <- stats::integrate(
      integrand,
      lower = 0,
      upper = upper,
      subdivisions = 2000,
      rel.tol = rel_tol,
      stop.on.error = FALSE
    )$value

    result <- 0.5 + integral_value

    min(max(result, 0), 1)
  }

  vapply(
    q,
    one_cdf,
    numeric(1)
  )
}



#' Quantiles of the Centered MHCP Log Process
#'
#' Computes quantiles of
#'
#' \deqn{
#' Z_H = \sum_{\ell=2}^{H} Y_\ell.
#' }
#'
#' @param p Probability in `(0,1)`.
#' @param H Component index. Must be at least 2.
#' @param rel_tol Numerical tolerance.
#'
#' @return Numeric quantile.
#'
#' @examples
#' qmhcp_centered(0.90, H = 10)
#'
#' @export
qmhcp_centered <- function(
    p,
    H,
    rel_tol = 1e-8
) {

  if (length(p) != 1L || p <= 0 || p >= 1) {
    stop("`p` must be a single probability strictly between 0 and 1.",
         call. = FALSE)
  }

  if (length(H) != 1L || H < 2 || H != as.integer(H)) {
    stop("`H` must be an integer >= 2.", call. = FALSE)
  }

  if (abs(p - 0.5) < 1e-14) {
    return(0)
  }

  # symmetry
  if (p < 0.5) {
    return(
      -qmhcp_centered(
        1 - p,
        H,
        rel_tol
      )
    )
  }

  sd_H <- (pi / 2) * sqrt(H - 1)

  upper <- 8 * sd_H

  while (
    pmhcp_centered(
      upper,
      H,
      rel_tol = rel_tol
    ) < p
  ) {
    upper <- 1.5 * upper
  }

  stats::uniroot(
    function(x) {
      pmhcp_centered(
        x,
        H,
        rel_tol = rel_tol
      ) - p
    },
    interval = c(0, upper),
    tol = rel_tol
  )$root
}





#' MHCP Marginal Distribution Function
#'
#' Computes `P(eta_H <= x)` for the MHCP.
#'
#' @param x Positive numeric value.
#' @param H Component index.
#' @param zeta Positive MHCP parameter.
#'
#' @return Numeric probability.
#'
#' @examples
#' pmhcp(
#'   x = 0.05,
#'   H = 10,
#'   zeta = 0.369411
#' )
#'
#' @export
pmhcp <- function(
    x,
    H,
    zeta
) {

  if (H < 1 || H != as.integer(H)) {
    stop("`H` must be a positive integer.", call. = FALSE)
  }

  if (zeta <= 0) {
    stop("`zeta` must be positive.", call. = FALSE)
  }

  result <- numeric(length(x))

  result[x <= 0] <- 0

  positive <- x > 0

  if (!any(positive)) {
    return(result)
  }

  if (H == 1) {

    result[positive] <- as.numeric(
      x[positive] >= 1
    )

    return(result)
  }

  threshold <- log(x[positive]) -
    (H - 1) * log(zeta)

  result[positive] <- pmhcp_centered(
    threshold,
    H = H
  )

  result
}
