
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mhcp

<!-- badges: start -->
<!-- badges: end -->

**`mhcp`** provides tools for simulation, probabilistic calibration, and
effective-dimension analysis of the **multiplicative half-Cauchy process
(MHCP)**.

The MHCP is an ordered shrinkage process defined by

$$
\eta_h = \prod_{\ell=1}^{h}\widetilde{\eta}_\ell,
$$

where

$$
\widetilde{\eta}_1 = 1,
\qquad
\widetilde{\eta}_h \sim C^+(0,\zeta),
\quad h \geq 2.
$$

On the log scale,

$$
\log \eta_h
=
(h-1)\log\zeta
+
\sum_{\ell=2}^{h}Y_\ell,
$$

where the centered increments $Y_\ell$ have density

$$
f_Y(y)=\frac{1}{\pi\cosh(y)}.
$$

This representation makes the MHCP useful as an ordered shrinkage prior:
when $0<\zeta<1$, trajectories shrink toward zero almost surely while
retaining the possibility of occasional large multiplicative excursions.

The package accompanies the manuscript:

> **Pathwise Shrinkage and Probabilistic Calibration of the
> Multiplicative Half-Cauchy Process**  
> Siddhesh Kulkarni (2026)

## Installation

You can install the development version of `mhcp` from GitHub with:

``` r
# install.packages("remotes")
remotes::install_github("SiddheshKulkarni-source/mhcp")
```

Then load the package:

``` r
library(mhcp)
```

## Quick start

### Simulate an MHCP

Use `rmhcp()` to generate trajectories from the process.

``` r
set.seed(123)

eta <- rmhcp(
  n = 5,
  d = 20,
  zeta = 0.7
)

round(eta[, 1:6], 3)
```

    ##      eta_1 eta_2 eta_3  eta_4  eta_5  eta_6
    ## [1,]     1 0.888 0.090  0.009  0.002  0.000
    ## [2,]     1 0.549 4.341 20.578 14.055 14.206
    ## [3,]     1 2.381 0.586  0.657  0.061  0.091
    ## [4,]     1 0.270 1.157  3.488  4.068  0.051
    ## [5,]     1 0.132 0.676  0.159  0.016  0.021

Each row represents an independent trajectory and each column
corresponds to an ordered component.

The first component satisfies $\eta_1=1$. Subsequent components are
generated multiplicatively.

For numerical work involving long trajectories, the process can also be
returned directly on the log scale:

``` r
log_eta <- rmhcp(
  n = 5,
  d = 20,
  zeta = 0.7,
  log = TRUE
)
```

## Pathwise shrinkage

A central property of the MHCP is

$$
\eta_h^{1/(h-1)}
\longrightarrow
\zeta
\qquad \text{almost surely}.
$$

Consequently,

$$
0<\zeta<1
\quad\Longrightarrow\quad
\eta_h\to0
\qquad \text{almost surely}.
$$

Thus, $\zeta$ determines the long-run geometric rate of the process.

Importantly, this does **not** mean that every realized trajectory is
monotonically decreasing. The half-Cauchy increments allow occasional
large upward movements.

## Local multiplicative excursions

The probability that the next component exceeds $c$ times the current
component has the closed form

$$
P(\eta_h>c\eta_{h-1})
=
\frac{2}{\pi}
\arctan\left(\frac{\zeta}{c}\right).
$$

The package evaluates this directly:

``` r
mhcp_excursion_prob(
  zeta = 0.7,
  c = c(1, 2, 5)
)
```

    ## [1] 0.38880022 0.21433385 0.08855123

For $\zeta=0.7$, this gives approximately

$$
P(\eta_h>\eta_{h-1})=0.389,
$$

while larger rebounds become progressively less likely.

This illustrates an important feature of the MHCP:

> **global shrinkage with local flexibility.**

## Probability-based prior calibration

One of the main purposes of `mhcp` is to translate an interpretable
prior statement into a value of the shrinkage parameter $\zeta$.

Suppose we want

$$
P(\eta_H\leq\varepsilon)=p.
$$

For example, suppose we want the scale of component 10 to be below
$0.05$ with 90% prior probability:

``` r
cal <- calibrate_mhcp(
  H = 10,
  epsilon = 0.05,
  p = 0.90
)

cal
```

    ##    H epsilon   p        q     zeta method
    ## 1 10    0.05 0.9 5.966875 0.369411  exact

The resulting value is approximately

$$
\zeta \approx 0.369.
$$

Thus,

$$
P(\eta_{10}\leq0.05)\approx0.90.
$$

We can verify this directly:

``` r
pmhcp(
  x = 0.05,
  H = 10,
  zeta = cal$zeta
)
```

    ## [1] 0.9

This should return a value close to

``` text
0.90
```

### Why this is useful

Instead of choosing $\zeta$ through trial and error, the analyst can
specify three interpretable quantities:

- $H$: the component at which shrinkage is desired,
- $\varepsilon$: a practically negligible scale,
- $p$: the desired prior probability of being below that scale.

The package then determines the corresponding $\zeta$.

The exact calibration is

$$
\zeta
=
\exp\left\{
\frac{\log\varepsilon-q_{p,H}}{H-1}
\right\},
$$

where $q_{p,H}$ is the $p$-quantile of the centered log process.

## Exact and approximate calibration

The centered log process

$$
Z_H=\sum_{\ell=2}^{H}Y_\ell
$$

has characteristic function

$$
\phi_{Z_H}(t)
=
\operatorname{sech}^{H-1}
\left(\frac{\pi t}{2}\right).
$$

`mhcp` evaluates its distribution using characteristic-function
inversion.

For example:

``` r
pmhcp_centered(
  q = 0,
  H = 10
)
```

    ## [1] 0.5

returns $0.5$, reflecting symmetry of the centered process.

Its quantiles can be obtained using:

``` r
qmhcp_centered(
  p = 0.90,
  H = 10
)
```

    ## [1] 5.966875

which gives approximately

$$
q_{0.90,10}\approx5.967.
$$

For faster approximate calculations, `calibrate_mhcp()` also provides a
normal approximation based on

$$
Z_H
\approx
N\left(
0,
\frac{(H-1)\pi^2}{4}
\right).
$$

``` r
calibrate_mhcp(
  H = 10,
  epsilon = 0.05,
  p = 0.90,
  method = "normal"
)
```

    ##    H epsilon   p        q      zeta method
    ## 1 10    0.05 0.9 6.039169 0.3664556 normal

## Threshold-based effective dimension

A calibration statement such as

$$
P(\eta_{10}\leq0.05)=0.90
$$

describes one selected component.

To summarize the entire sequence, define the threshold-based effective
dimension

$$
N_\varepsilon
=
\#\{h:\eta_h>\varepsilon\}.
$$

For $0<\zeta<1$,

$$
N_\varepsilon<\infty
\qquad \text{almost surely}.
$$

The package numerically evaluates its expectation:

``` r
eff <- expected_mhcp_dimension(
  zeta = cal$zeta,
  epsilon = 0.05
)

eff$mean_effective_dimension
```

    ## [1] 4.752268

For the calibration

$$
H=10,\qquad
\varepsilon=0.05,\qquad
p=0.90,
$$

the expected threshold-based effective dimension is approximately

$$
E(N_{0.05})\approx4.75.
$$

This provides a useful bridge between a local calibration statement and
the global complexity implied by the prior.

## Effect of the calibration probability

The probability level $p$ controls how aggressively the process is
shrunk.

``` r
p_values <- c(0.50, 0.75, 0.90)

calibration_table <- do.call(
  rbind,
  lapply(
    p_values,
    function(p) {
      calibrate_mhcp(
        H = 10,
        epsilon = 0.05,
        p = p
      )
    }
  )
)

calibration_table
```

    ##    H epsilon    p        q      zeta method
    ## 1 10    0.05 0.50 0.000000 0.7168712  exact
    ## 2 10    0.05 0.75 3.105284 0.5076867  exact
    ## 3 10    0.05 0.90 5.966875 0.3694110  exact

For $H=10$ and $\varepsilon=0.05$, the calibrated values are
approximately:

|  $p$ | $\zeta$ | $E(N_{0.05})$ |
|-----:|--------:|--------------:|
| 0.50 |   0.717 |         20.63 |
| 0.75 |   0.508 |          7.60 |
| 0.90 |   0.369 |          4.75 |

Higher values of $p$ correspond to stronger prior confidence that the
target component is negligible and therefore imply more aggressive
shrinkage.

## Visualizing MHCP trajectories

The following example illustrates the pathwise behavior of the process.

``` r
set.seed(123)

sim <- rmhcp(
  n = 30,
  d = 30,
  zeta = 0.7,
  log = TRUE
)

matplot(
  x = seq_len(ncol(sim)),
  y = t(sim),
  type = "l",
  lty = 1,
  xlab = "Component index h",
  ylab = expression(log(eta[h]))
)

lines(
  seq_len(ncol(sim)),
  (seq_len(ncol(sim)) - 1) * log(0.7),
  lwd = 3,
  lty = 2
)

legend(
  "topright",
  legend = c(
    "Sample trajectories",
    "Theoretical median"
  ),
  lty = c(1, 2),
  lwd = c(1, 3),
  bty = "n"
)
```

![](README_files/figure-gfm/trajectories-1.png)<!-- -->

Individual trajectories can exhibit substantial local variation, but the
long-run log-scale trend is governed by

$$
(h-1)\log\zeta.
$$

## Main functions

The current package interface is intentionally small.

| Function                    | Purpose                                              |
|:----------------------------|:-----------------------------------------------------|
| `rmhcp()`                   | Simulate MHCP trajectories                           |
| `pmhcp()`                   | Evaluate $P(\eta_H\leq x)$                           |
| `pmhcp_centered()`          | Evaluate the CDF of the centered log process         |
| `qmhcp_centered()`          | Compute centered log-process quantiles               |
| `calibrate_mhcp()`          | Calibrate $\zeta$ from $P(\eta_H\leq\varepsilon)=p$  |
| `mhcp_excursion_prob()`     | Compute local multiplicative excursion probabilities |
| `expected_mhcp_dimension()` | Compute $E(N_\varepsilon)$                           |

Individual help pages are available through R, for example:

``` r
?calibrate_mhcp
?rmhcp
?expected_mhcp_dimension
```

## Reproducibility

The numerical tools in this package are designed to reproduce the
calibration and stochastic-process calculations described in the
accompanying manuscript.

To reproduce the main worked calibration:

``` r
result <- calibrate_mhcp(
  H = 10,
  epsilon = 0.05,
  p = 0.90
)

result
```

    ##    H epsilon   p        q     zeta method
    ## 1 10    0.05 0.9 5.966875 0.369411  exact

``` r
expected_mhcp_dimension(
  zeta = result$zeta,
  epsilon = 0.05
)$mean_effective_dimension
```

    ## [1] 4.752268

## Citation

If you use `mhcp` in research, please cite the accompanying manuscript:

> Kulkarni, S. (2026). *Pathwise Shrinkage and Probabilistic Calibration
> of the Multiplicative Half-Cauchy Process*. Manuscript in preparation.

The package citation can also be displayed in R using:

``` r
citation("mhcp")
```

The citation will be updated when the accompanying article receives its
final bibliographic information.

## Development status

`mhcp` is currently under active development.

The current version focuses on the theoretical and computational tools
needed for simulation, probabilistic prior calibration, excursion
probabilities, and threshold-based effective-dimension analysis.

Bug reports and suggestions are welcome through the GitHub issue
tracker.

## License

`mhcp` is released under the MIT License.
