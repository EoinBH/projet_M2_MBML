# coclustOrdinal

Co-clustering of ordinal data (for example questionnaire answers on a
1-to-*m* scale) with a **latent block model** whose blocks follow a **CUB**
distribution. Parameters are estimated with a **SEM-Gibbs** algorithm and the
numbers of row and column clusters are chosen with the **ICL** criterion.

*Projet du cours Model-Based Machine Learning, M2 MALIA, Université Lumière
Lyon 2 (2026-2027).*

## Installation

The package is written in pure R (no compilation needed) and only depends on
`stats`. From the source tarball:

```r
install.packages("coclustOrdinal_0.1.0.tar.gz", repos = NULL, type = "source")
```

or from a clone of this repository:

```r
# install.packages("devtools")
devtools::install(build_vignettes = TRUE)
```

Building the vignette needs `knitr` and `rmarkdown`. It also uses `mclust`
(optional) for the adjusted Rand index.

## Quick start

```r
library(coclustOrdinal)

# Simulated 100 x 50 matrix with 2 row clusters and 2 column clusters
sim <- simulate_cub_lbm(n = 100, d = 50, m = 5, seed = 1)

# Fit with K = 2 row clusters and L = 2 column clusters (20 random starts)
fit <- coclust_ordinal(sim$X, K = 2, L = 2, seed = 1)
summary(fit)

fit$row_clusters    # row assignments
fit$col_prob        # column membership probabilities
fit$icl             # ICL of the fit
fit$params          # alpha, beta, pi, xi

# Choose K and L by ICL
sel <- coclust_ordinal_select(sim$X, K_range = 1:4, L_range = 1:4, n_init = 5)
sel$table
```

## Functions

| Function | Purpose |
|---|---|
| `coclust_ordinal()` | Fit the CUB latent block model for given `K` and `L`, from `n_init` random initializations (20 by default), and return the best solution by ICL. |
| `coclust_ordinal_select()` | Fit every `(K, L)` of a grid and select the pair with the largest ICL. |
| `simulate_cub_lbm()` | Simulate ordinal data from the model, with the true partitions. |

`coclust_ordinal()` returns the row and column membership probabilities
(`row_prob`, `col_prob`), the assignments (`row_clusters`, `col_clusters`),
the ICL value (`icl`) and the estimated parameters (`params`).

## Documentation

```r
?coclust_ordinal
browseVignettes("coclustOrdinal")
```

The vignette explains the model, the algorithm and the ICL criterion, with a
worked example.

## Code organisation

| File | Contents |
|---|---|
| `R/coclust_ordinal.R` | Main functions `coclust_ordinal()` and `coclust_ordinal_select()` |
| `R/sem_gibbs.R` | SEM-Gibbs algorithm: main loop, initialization, Gibbs sweep, M step |
| `R/cub.R` | CUB probabilities and EM estimation of `(pi, xi)` in one block |
| `R/icl.R` | Complete-data log-likelihood and ICL |
| `R/utils.R` | Input checks and small matrix helpers |
| `R/simulate.R` | Data simulation |
| `R/methods.R` | `print()` and `summary()` methods |

## References

- D. Piccolo (2003). On the moments of a mixture of uniform and shifted
  binomial random variables. *Quaderni di Statistica*, 5. (The CUB model.)
- J. Jacques and C. Biernacki (2018). Model-based co-clustering for ordinal
  data. *Computational Statistics & Data Analysis*. (Latent block model,
  SEM-Gibbs and ICL for ordinal data, there with the BOS distribution
  instead of CUB.)

## Author

Eoin Brereton Hurley
