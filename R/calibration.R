#' Calibrate the MHCP Shrinkage Parameter
#'
#' Chooses `zeta` so that
#'
#' \deqn{
#' P(\eta_H \le \epsilon) = p.
#' }
#'
#' The semi-analytic calibration is
#'
#' \deqn{
#' \zeta
#' =
#' \exp\left[
#' \frac{
#' \log(\epsilon)-q_{p,H}
#' }{
#' H-1
#' }
#' \right],
#' }
#'
#' where `q_{p,H}` is the corresponding quantile of the centered
#' log process.
#'
#' @param H Target component index.
#' @param epsilon Positive shrinkage threshold.
#' @param p Target probability.
#' @param method Either `"exact"` for characteristic-function inversion
#'   or `"normal"` for the CLT approximation.
#'
#' @return A data frame containing the calibration inputs and `zeta`.
#'
#' @examples
#' calibrate_mhcp(
#'   H = 10,
#'   epsilon = 0.05,
#'   p = 0.90
#' )
#'
#' @export
calibrate_mhcp <- function(
    H,
    epsilon,
    p,
    method = c("exact", "normal")
) {

  method <- match.arg(method)

  if (H < 2 || H != as.integer(H)) {
    stop("`H` must be an integer >= 2.", call. = FALSE)
  }

  if (!is.finite(epsilon) || epsilon <= 0) {
    stop("`epsilon` must be positive.", call. = FALSE)
  }

  if (!is.finite(p) || p <= 0 || p >= 1) {
    stop("`p` must lie strictly between 0 and 1.", call. = FALSE)
  }

  if (method == "exact") {

    q_value <- qmhcp_centered(
      p = p,
      H = H
    )

  } else {

    q_value <- (pi / 2) *
      sqrt(H - 1) *
      stats::qnorm(p)
  }

  zeta <- exp(
    (
      log(epsilon) -
        q_value
    ) /
      (H - 1)
  )

  data.frame(
    H = H,
    epsilon = epsilon,
    p = p,
    q = q_value,
    zeta = zeta,
    method = method
  )
}
