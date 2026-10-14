# MC-PLS for multilevel SEMs

Simulation study comparing the multilevel MC-PLS estimator in `plssem`
(`pls(..., cluster = )`) with Mplus BAYES (via `modsem::modsem_mplus()`),
using continuous and ordinal indicators.

## Population models

Two two-level SEMs, based on Cox & Kelcey (2025). See also Cox, Kelcey & Bai
(2023). Both models share the same measurement model and level-2 model.

The loadings and ICCs are the standardized equivalents of Cox & Kelcey's
model (their Table 1). In their model, the X and Z factors have unit variance,
the loadings are 1, the error variances are .8 within and .2 between, and Y
has residual variances of .8 within and .2 between.

| Parameter                    | Within (level 1) | Between (level 2)    |
|------------------------------|------------------|----------------------|
| Loadings X, Z (3 per factor) | .745 (ω = .79)   | .913 (ω = .94)       |
| Loadings Y (3)               | .780             | .866                 |
| Covariance                   | `Xw ~~ .25 Zw`   | `Xb ~~ .25 Zb`       |
| Structural                   | (see below)      | `Yb ~ .4 Xb + .4 Zb` |

The models differ in the within-level structural model:

1. **`within_xz`** (within-level latent interaction):
   `Yw ~ .4 Xw + .4 Zw + .2 Xw:Zw`.
2. **`random_slope`** (cross-level interaction through a random slope):
   `Yw ~ .4 Xw + rv("s") Xw + .4 Zw`, with `s ~ .5 Zb` at level 2.
   - The slope SD is `rsd = .4`, so the slope in cluster j is `.4 + .4 s_j`.
   - The implied cross-level interaction is .4 × .5 = .2, as in Cox & Kelcey.
   - The random slope is written as a duplicated statement. The fixed effect
     is removed from the estimation syntax by `get_estimation_syntax()`.

- The ICC is .40 for the x and z indicators and .281 for the y indicators.
- All parameters are standardized at each level. The structural model is thus
  not identical to Cox & Kelcey's: there, Y doesn't have unit variance, so
  their effects of .4 and .2 correspond to .36 and .18 (within) and .52
  (between) on the standardized scale.
- The latent variables have different names at each level (`*w`/`*b`), as
  Mplus requires.

## Design

| Factor                 | Levels                                          |
|------------------------|-------------------------------------------------|
| Model                  | `within_xz`, `random_slope`                     |
| Cluster size (n)       | 20, 50                                          |
| Number of clusters (k) | 30, 60, 90                                      |
| Data                   | continuous                                      |
|                        | ordinal: 2 or 5 categories × Symmetric/Moderate |

- This gives 2 × 6 × 5 = 60 cells, with R = 200 replications per `run.id`.
- n and k follow Cox, Kelcey & Bai (2023) and Cox & Kelcey (2025).
- The thresholds follow Rhemtulla et al. (2012) and Schubert et al. (2018).

## Estimators

- **MC-PLS:** `plssem::pls()` with 500 bootstrap samples for the SEs.
- **Mplus:** `ESTIMATOR = BAYES` with default priors and
  `BITERATIONS = (20000)`, as in Cox & Kelcey (2025).

Both standardized and unstandardized estimates are recorded. True values are
only defined for the standardized ones, so `true` is `NA` for the
unstandardized rows. Level-2 parameters have the suffix `.l2`.

The slope SD (`rsd.s`) is only recorded for MC-PLS, because Mplus only
reports the (unstandardized) slope variance. Mplus reports the fixed effect of
the random slope as `s ~1`, so the unstandardized `Yw~Xw` is `NA` for Mplus.
The standardized `Yw~Xw` is reported by both estimators.

## References

- Cox, K., & Kelcey, B. (2025). Comparing Bayesian estimation and
  structural-after-measurement approaches for structural equation models with
  latent interactions and complex data structures. *Behavior Research Methods,
  57*, 320. https://doi.org/10.3758/s13428-025-02840-x
- Cox, K., Kelcey, B., & Bai, ?. (2023). Croon's bias-corrected estimation for
  multilevel structural equation models with latent interactions. *Structural
  Equation Modeling, 30*(3), 467–480.
- Rhemtulla, M., Brosseau-Liard, P. É., & Savalei, V. (2012). When can
  categorical variables be treated as continuous? *Psychological Methods,
  17*(3), 354–373. https://doi.org/10.1037/a0029315
