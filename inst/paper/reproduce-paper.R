library(mhcp)

## ============================================================
## MHCP paper reproducibility script
## ============================================================

H <- 10
epsilon <- 0.05
p_values <- c(0.50, 0.75, 0.90)

## ------------------------------------------------------------
## 1. Calibration results
## ------------------------------------------------------------

results <- do.call(
  rbind,
  lapply(
    p_values,
    function(p) {

      cal <- calibrate_mhcp(
        H = H,
        epsilon = epsilon,
        p = p
      )

      eff <- expected_mhcp_dimension(
        zeta = cal$zeta,
        epsilon = epsilon
      )

      data.frame(
        H = H,
        epsilon = epsilon,
        p = p,
        q = cal$q,
        zeta = cal$zeta,
        expected_dimension = eff$mean_effective_dimension
      )
    }
  )
)

print(results, digits = 8)

## ------------------------------------------------------------
## 2. Verify calibration numerically
## ------------------------------------------------------------

results$verified_probability <- mapply(
  function(zeta, p) {
    pmhcp(
      x = epsilon,
      H = H,
      zeta = zeta
    )
  },
  zeta = results$zeta,
  p = results$p
)

print(results, digits = 8)

## ------------------------------------------------------------
## 3. Local excursion probabilities
## ------------------------------------------------------------

excursions <- data.frame(
  c = c(1, 2, 5),
  probability = mhcp_excursion_prob(
    zeta = 0.7,
    c = c(1, 2, 5)
  )
)

print(excursions, digits = 8)

## ------------------------------------------------------------
## 4. Median calibration identity
## ------------------------------------------------------------

median_exact <- calibrate_mhcp(
  H = H,
  epsilon = epsilon,
  p = 0.50
)$zeta

median_closed_form <- epsilon^(1 / (H - 1))

cat("\nMedian calibration\n")
cat("Package:     ", median_exact, "\n")
cat("Closed form: ", median_closed_form, "\n")

stopifnot(
  isTRUE(
    all.equal(
      median_exact,
      median_closed_form,
      tolerance = 1e-7
    )
  )
)

## ------------------------------------------------------------
## 5. Figure: pathwise behavior
## ------------------------------------------------------------

set.seed(2026)

zeta_values <- c(0.7, 1.0, 1.3)
d <- 40
n_paths <- 20

old_par <- par(
  mfrow = c(1, 3),
  mar = c(4, 4, 3, 1)
)

for (zeta in zeta_values) {

  sim <- rmhcp(
    n = n_paths,
    d = d,
    zeta = zeta,
    log = TRUE
  )

  h <- seq_len(d)

  matplot(
    h,
    t(sim),
    type = "l",
    lty = 1,
    xlab = "Component index h",
    ylab = "log(eta_h)",
    main = paste("zeta =", zeta)
  )

  lines(
    h,
    (h - 1) * log(zeta),
    lwd = 3,
    lty = 2
  )
}

par(old_par)

## ------------------------------------------------------------
## 6. Session information
## ------------------------------------------------------------

sessionInfo()
