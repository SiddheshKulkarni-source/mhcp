# Internal numerically stable hyperbolic secant
sech <- function(x) {

  ax <- abs(x)

  2 * exp(-ax) /
    (1 + exp(-2 * ax))
}
