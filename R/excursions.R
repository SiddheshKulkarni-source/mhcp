#' MHCP Multiplicative Excursion Probability
#'
#' Computes
#'
#' \deqn{
#' P(\eta_h > c\eta_{h-1})
#' =
#' \frac{2}{\pi}
#' \arctan\left(\frac{\zeta}{c}\right).
#' }
#'
#' @param zeta Positive MHCP parameter.
#' @param c Positive multiplicative excursion size.
#'
#' @return Numeric probability.
#'
#' @examples
#' mhcp_excursion_prob(0.7)
#' mhcp_excursion_prob(0.7, c = 5)
#'
#' @export
mhcp_excursion_prob <- function(
    zeta,
    c = 1
) {

  if (any(zeta <= 0)) {
    stop("`zeta` must be positive.", call. = FALSE)
  }

  if (any(c <= 0)) {
    stop("`c` must be positive.", call. = FALSE)
  }

  (2 / pi) *
    atan(
      zeta / c
    )
}
