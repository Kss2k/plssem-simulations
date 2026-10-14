# ──────────────────────────────────────────────────────────────────────────────
# Setup
# ──────────────────────────────────────────────────────────────────────────────
devtools::load_all() # load custom utils

library(mvtnorm)
library(tidyr)
library(dplyr)
library(modsem) # v.1.0.24 (main, >= cbdb438b, needs multilevel `modsem_mplus()`)
library(plssem) # v.0.2.0

.root_info   <- set_project_root()
setwd("mcpls-mlm")

# Absolute path, so parallel workers (which do not inherit the master's working
# directory) read/write results in the right place regardless of their wd.
PROJECT_ROOT <- .root_info$project.root
RESULTS_DIR <- normalizePath("results", mustWork = FALSE)
dir.create(RESULTS_DIR, showWarnings = FALSE)

# Log file for this session. A new file is created at each startup (so the logs
# of previous sessions are kept), and each (parallel) worker appends a line for
# every finished iteration, as the live output of the workers isn't visible when
# running in parallel. Follow it with e.g. `tail -f logs/log-<stamp>.txt`.
LOG_DIR  <- normalizePath("logs", mustWork = FALSE)
dir.create(LOG_DIR, showWarnings = FALSE)
LOG_FILE <- file.path(LOG_DIR, sprintf("log-%s.txt", format(Sys.time(), "%Y-%m-%d-%H-%M-%S")))

write_log <- function(...) {
  line <- sprintf("[%s | pid %d] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                  Sys.getpid(), paste0(...))
  cat(line, "\n", file = LOG_FILE, append = TRUE, sep = "")
}

write_log("Session started")

# Run settings. Edit these before sourcing/running this script in R.
checkIfExists <- TRUE
R             <- 200L
run.id        <- NULL

parallel  <- TRUE
n.workers <- 8

# ──────────────────────────────────────────────────────────────────────────────
# Model+Parameters
# ──────────────────────────────────────────────────────────────────────────────

# Measurement model and ICCs from Cox & Kelcey (2025, Table 1). In their
# (unstandardized) model, Xw, Xb, Zw and Zb have unit variances, the indicators
# have loadings of 1.0 and error variances of 0.8 (within) and 0.2 (between),
# and Y has residual variances of 0.8 (within) and 0.2 (between). With
# beta = 0.4, beta_XZ = 0.2 and cov(X, Z) = 0.25, the variance of Y is
#   within:  0.4^2 + 0.4^2 + 2 * 0.4^2 * 0.25 + 0.2^2 * (1 + 0.25^2) + 0.8 = 1.2425
#   between: 0.4^2 + 0.4^2 + 2 * 0.4^2 * 0.25 + 0.2                    = 0.6
# The standardized loadings (each level is standardized here) are thus:
1 / sqrt(1 + 0.8)                # X, Z within
#> [1] 0.745356
1 / sqrt(1 + 0.2)                # X, Z between
#> [1] 0.9128709
sqrt(1.2425 / (1.2425 + 0.8))    # Y within
#> [1] 0.7799507
sqrt(0.6 / (0.6 + 0.2))          # Y between
#> [1] 0.8660254
# and the ICCs of the indicators are:
1.2 / (1.8 + 1.2)                # x, z
#> [1] 0.4
0.8 / ((1.2425 + 0.8) + 0.8)     # y
#> [1] 0.2814424

# Formula for composite reliability:
cr <- function(...) (sum(c(...))^2)/(sum(c(...))^2 + sum(1-c(...)^2))
cr(0.745, 0.745, 0.745) # X, Z within (omega = 0.79 in Cox & Kelcey)
#> [1] 0.7891164
cr(0.913, 0.913, 0.913) # X, Z between (omega = 0.94 in Cox & Kelcey)
#> [1] 0.9375994

icc <- c(
  x1 = 0.400, x2 = 0.400, x3 = 0.400,
  z1 = 0.400, z2 = 0.400, z3 = 0.400,
  y1 = 0.281, y2 = 0.281, y3 = 0.281
)

# Model 1 (within_xz):
# Based on the within-level latent interaction model in Cox & Kelcey (2025),
# see also Cox, Kelcey & Bai (2023). The measurement model and ICCs are the
# same as in Cox & Kelcey (see above). The structural model is standardized at
# each level (unlike in Cox & Kelcey, where Y doesn't have a unit variance):
# path coefficients are 0.4, the interaction is 0.2, and the covariance between
# X and Z is 0.25 at both levels.
#
# Model 2 (random_slope):
# Same measurement model and level 2 model as model 1, but instead of the
# within-level interaction the effect of Xw on Yw is a random slope `s`, which
# is predicted by Zb (i.e., a cross-level interaction, similar to the
# cross-level model in Cox & Kelcey, 2025, but with residual slope variance).
# The fixed effect is 0.4, and the random slope is specified with a duplicated
# statement (`rv("s") * Xw`, see `get_estimation_syntax()`). `s` is
# standardized at level 2, and the slope in cluster j is `0.4 + rsd * s_j`,
# with `rsd = 0.4` and `s ~ 0.5 * Zb`. The implied cross-level interaction is
# thus 0.4 * 0.5 = 0.2, as in Cox & Kelcey (2025), and Zb explains 25% of the
# variance of the slope.
#
# Mplus requires different names for the latent variables at each level, so
# we use `*w` (within) and `*b` (between).
models <- list(
  within_xz = list(
    syntax = '
      level: 1
         Xw =~ 0.745 * x1 + 0.745 * x2 + 0.745 * x3
         Zw =~ 0.745 * z1 + 0.745 * z2 + 0.745 * z3
         Yw =~ 0.780 * y1 + 0.780 * y2 + 0.780 * y3
         Yw  ~ 0.4 * Xw + 0.4 * Zw + 0.2 * Xw:Zw
         Xw ~~ 0.25 * Zw
      level: 2
         Xb =~ 0.913 * x1 + 0.913 * x2 + 0.913 * x3
         Zb =~ 0.913 * z1 + 0.913 * z2 + 0.913 * z3
         Yb =~ 0.866 * y1 + 0.866 * y2 + 0.866 * y3
         Yb  ~ 0.4 * Xb + 0.4 * Zb
         Xb ~~ 0.25 * Zb
    ',
    icc = icc
  ),

  random_slope = list(
    syntax = '
      level: 1
         Xw =~ 0.745 * x1 + 0.745 * x2 + 0.745 * x3
         Zw =~ 0.745 * z1 + 0.745 * z2 + 0.745 * z3
         Yw =~ 0.780 * y1 + 0.780 * y2 + 0.780 * y3
         Yw  ~ 0.4 * Xw + rv("s") * Xw + 0.4 * Zw
         Xw ~~ 0.25 * Zw
      level: 2
         Xb =~ 0.913 * x1 + 0.913 * x2 + 0.913 * x3
         Zb =~ 0.913 * z1 + 0.913 * z2 + 0.913 * z3
         Yb =~ 0.866 * y1 + 0.866 * y2 + 0.866 * y3
         Yb  ~ 0.4 * Xb + 0.4 * Zb
         s   ~ 0.5 * Zb
         Xb ~~ 0.25 * Zb
    ',
    icc = icc,
    rsd = c(s = 0.4)
  )
)

# Estimation syntax, and the true values of each model
for (m in names(models)) {
  models[[m]]$syntax.est <- get_estimation_syntax(models[[m]]$syntax)
}


# The true values are only known for the standardized parameters, but
# we also record the unstandardized estimates. Parameters at level 2
# get a `.l2` suffix (as in `plssem`). The standard deviations of the random
# slopes (`rsd`) are only estimated by MC-PLS (`rsd.<name>`), and are only
# recorded as standardized parameters.
get_true_values <- function(syntax, rsd = NULL) {
  levels <- split_levels(modsemify(syntax))

  std <- do.call(rbind, lapply(1:2, FUN = \(l) {
    pt <- levels[[l]]
    pt <- pt[pt$op %in% c("=~", "~", "~~") & !is_random_slope(pt$mod), , drop = FALSE]

    data.frame(
      par   = paste0(pt$lhs, pt$op, pt$rhs, if (l == 2L) ".l2" else ""),
      true  = as.numeric(pt$mod),
      level = l
    )
  }))

  ustd <- std
  ustd$true <- NA_real_

  out <- rbind(
    cbind(std, type = "standardized"),
    cbind(ustd, type = "unstandardized")
  )

  if (length(rsd)) {
    out <- rbind(out, data.frame(
      par   = paste0("rsd.", names(rsd)),
      true  = unname(rsd),
      level = 2L,
      type  = "standardized"
    ))
  }

  out$key <- paste(out$type, out$par, sep = ":")
  out
}


for (m in names(models)) {
  models[[m]]$true <- get_true_values(models[[m]]$syntax, rsd = models[[m]]$rsd)
}


# Based on Rhemtulla et al., 2012 and Schubert et al., 2018
list_thresholds <- list(
  Symmetric = list(
    `2` = c( 0.00),
    `3` = c(-0.83,  0.83),
    `4` = c(-1.25,  0.00,  1.25),
    `5` = c(-1.50, -0.50,  0.50, 1.50),
    `6` = c(-1.60, -0.83,  0.00, 0.83, 1.60),
    `7` = c(-1.79, -1.07, -0.36, 0.36, 1.07, 1.79)
  ),
  Moderate = list(
    `2` = c( 0.36),
    `3` = c(-0.50,  0.76),
    `4` = c(-0.31,  0.79,  1.66),
    `5` = c(-0.70,  0.39,  1.16,  2.05),
    `6` = c(-1.05,  0.08,  0.81,  1.44,  2.33),
    `7` = c(-1.43, -0.43,  0.38,  0.94,  1.44,  2.54)
  ),
  Extreme = list(
    `2` = c( 1.04),
    `3` = c( 0.58,  1.13),
    `4` = c( 0.28,  0.71,  1.23),
    `5` = c( 0.05,  0.44,  0.84,  1.34),
    `6` = c(-0.13,  0.25,  0.61,  0.99,  1.48),
    `7` = c(-0.25,  0.13,  0.47,  0.81,  1.18,  1.64)
  ),
  Alt.Mod = list(
    `2` = c(-0.36),
    `3` = c(-0.76,  0.50),
    `4` = c(-1.66, -0.79,  0.31),
    `5` = c(-2.05, -1.16, -0.39,  0.70),
    `6` = c(-2.33, -1.44, -0.81, -0.08,  1.05),
    `7` = c(-2.54, -1.44, -0.94, -0.38,  0.43,  1.43)
  ),
  Alt.Ext = list(
    `2` = c(-1.04),
    `3` = c(-1.13, -0.58),
    `4` = c(-1.23, -0.71, -0.28),
    `5` = c(-1.34, -0.84, -0.44, -0.05),
    `6` = c(-1.48, -0.99, -0.61, -0.25,  0.13),
    `7` = c(-1.64, -1.18, -0.81, -0.47, -0.13,  0.25)
  )
)


# sample size based on Cox, Kelcey & Bai (2023)
n <- c(20, 50) # within cluster sample size
k <- c(30, 60, 90) # cluster sample size

# Data conditions. Continuous data (ncat = "Inf") is not crossed with skew.
# The ordinal conditions cover binary, and five-point (Likert type) data, with
# and without skew.
data.conditions <- data.frame(
  ncat = c("Inf",  "2",         "2",        "5",         "5"),
  skew = c("None", "Symmetric", "Moderate", "Symmetric", "Moderate")
)

# Set up selection indices which are crossed
idx.model <- seq_along(models)
idx.n     <- seq_along(n)
idx.k     <- seq_along(k)
idx.data  <- seq_len(NROW(data.conditions))

IDX <- expand.grid(
  model = idx.model,
  n     = idx.n,
  k     = idx.k,
  data  = idx.data
)

IDX$ncat <- data.conditions$ncat[IDX$data]
IDX$skew <- data.conditions$skew[IDX$data]

# ──────────────────────────────────────────────────────────────────────────────
# Estimators
# ──────────────────────────────────────────────────────────────────────────────

# Keys for matching the estimates to the true values, see `get_true_values()`
par_keys <- function(par, type) {
  paste0(type, ":", par$lhs, par$op, par$rhs, ifelse(par$level %in% 2, ".l2", ""))
}


est_pls <- function(model, data, ...) {
  fit <- plssem::pls(model, data, ...)

  spar <- as.data.frame(plssem::parameter_estimates(fit))
  upar <- as.data.frame(plssem::unstandardized_estimates(fit))

  coef <- rbind(
    cbind(spar$est, spar$se),
    cbind(upar$est, upar$se)
  )

  rownames(coef) <- c(par_keys(spar, "standardized"), par_keys(upar, "unstandardized"))
  colnames(coef) <- c("est", "se")

  # standard deviations of the random slopes (if any)
  rsd <- fit@params$rsd
  if (length(rsd)) {
    rsd.se <- fit@params$rsd.se
    if (is.null(rsd.se)) rsd.se <- rep(NA_real_, length(rsd))

    rsd <- cbind(est = unname(rsd), se = unname(rsd.se[names(rsd)]))
    rownames(rsd) <- paste0("standardized:rsd.", names(fit@params$rsd))
    coef <- rbind(coef, rsd)
  }

  # NOTE: `checkIfParTableIsAdmissible()` doesn't work for multilevel models
  attr(coef, "admissible") <- is_admissible(fit)

  coef
}


est_mplus <- function(model, data, ...) {
  fit <- modsem::modsem_mplus(model, data, ...)

  upar <- as.data.frame(modsem::parameter_estimates(fit))
  spar <- as.data.frame(modsem::standardized_estimates(fit))

  # Mplus is case insensitive, so we map the names back to the
  # (case sensitive) names in the model
  mod   <- modsemify(model)
  vars0 <- setdiff(union(mod$lhs, mod$rhs), "level")
  remap <- function(x) {
    y <- vars0[match(tolower(x), tolower(vars0))]
    ifelse(is.na(y), x, y)
  }

  for (nm in c("lhs", "rhs")) {
    upar[[nm]] <- remap(upar[[nm]])
    spar[[nm]] <- remap(spar[[nm]])
  }

  coef <- rbind(
    cbind(spar$est, spar$std.error),
    cbind(upar$est, upar$std.error)
  )

  rownames(coef) <- c(par_keys(spar, "standardized"), par_keys(upar, "unstandardized"))
  colnames(coef) <- c("est", "se")

  # Non-converged models don't produce estimates (i.e., `modsem_mplus()` fails)
  attr(coef, "admissible") <- all(is.finite(spar$est))

  coef
}


run_estimators <- function(methods,
                           data_i,
                           ordered,
                           model,
                           id,
                           n.i,
                           k.i,
                           skew,
                           ncat,
                           model.id,
                           seed) {
  methods <- match.arg(
    methods,
    choices = c("MC-PLS", "Mplus"),
    several.ok = TRUE
  )

  run_one <- function(method) {
    switch(
      method,
      "MC-PLS" = get_output(
        func          = est_pls,
        data          = data_i,
        cluster       = "cluster",
        ordered       = ordered,
        bootstrap     = TRUE,
        boot.R        = 500,
        model         = model$syntax.est,
        true          = model$true,
        extra         = list(k = k.i),
        method        = "MC-PLS",
        boot.ncores   = 2,
        boot.parallel = "multicore",
        verbose       = FALSE,
        id            = id,
        n             = n.i,
        skew          = skew,
        ncat          = ncat,
        model.id      = model.id,
        seed          = seed
      ),

      "Mplus" = get_output(
        func             = est_mplus,
        data             = data_i,
        cluster          = "cluster",
        model            = model$syntax.est,
        true             = model$true,
        extra            = list(k = k.i),
        cleanup          = TRUE,
        method           = "Mplus",
        processors       = 2,
        categorical      = ordered,
        estimator        = "bayes",
        analysis.options = list(BITERATIONS = "(20000)"),
        n                = n.i,
        id               = id,
        skew             = skew,
        ncat             = ncat,
        model.id         = model.id,
        seed             = seed
      )
    )
  }

  stats::setNames(lapply(methods, run_one), methods)
}


drop_csv_row_names <- function(x) {
  if (NROW(x) && names(x)[[1L]] %in% c("", "X", "X.1")) {
    x <- x[-1L]
  }

  x
}


# ──────────────────────────────────────────────────────────────────────────────
# Run Simulation
# ──────────────────────────────────────────────────────────────────────────────
K             <- NROW(IDX)
total         <- R * K


# ── Parallelisation ───────────────────────────────────────────────────────────
# Set `parallel <- TRUE` to evaluate the R batches concurrently with the
# `future` package. Each batch (one value of the outer `i` index) is fully
# independent: it derives its iteration ids from `i`, sets its own per-iteration
# seeds, and writes its own results CSV. The batches can therefore run in any
# order / in parallel and still reproduce the sequential output exactly.

# The run.id specifies the circumstance the script is running under
# v0-test is for testing. The other run.ids (see below) specify
# different computers running parallel simulations. They of course
# need different seeds, to generate unique results.

LOCAL_SEEDS <- c(
  "v1-test"    = 5340956,
  "v1-tuf"     = 1210967,
  "v1-promax"  = 2983429,
  "v1-rebekka" = 5700825
)

if (is.null(run.id)) {
  cat("What run.id do you want to use? Available:\n")
  print(names(LOCAL_SEEDS))
  run.id.idx <- as.integer(readLines(n=1))
  run.id <- names(LOCAL_SEEDS)[[run.id.idx]]
}

write_log(sprintf(
  "run.id = %s, R = %d, cells = %d, iterations = %d, parallel = %s, workers = %d",
  run.id, R, K, total, parallel, if (parallel) n.workers else 1L
))

# The run.id specifies what seed we set
set.seed(LOCAL_SEEDS[[run.id]])

# each iteration has it own seed, such that we can reproduce a specific
# iterartion in isolation (if desired). This seed is appended to the output.
seeds <- floor(runif(total, min = 0, max = 9999999))

# Run a single batch `i` (one full pass over IDX) and return its results.
# Self-contained so it can be called sequentially or dispatched to a future
# worker: ids and seeds are derived from `i`, and it writes its own CSV.
run_batch <- function(i) {
  results.i  <- NULL
  filePrefix <- paste("results", run.id, i, sep = "-")

  files <- dir(RESULTS_DIR)
  match <- startsWith(files, paste0(filePrefix, "-"))

  if (checkIfExists && any(match)) {
    message(sprintf(
      "Skipping iteration batch %d, as it has already been run...", i)
    )
    write_log(sprintf("batch %d: skipped (already run)", i))

    results.i <- read.csv(file.path(RESULTS_DIR, files[which(match)[[1L]]]))
    results.i <- results.i[-1] # drop rownames
    return(results.i)
  }

  for (j in seq_len(NROW(IDX))) {

    idx.modj  <- IDX$model[[j]]
    idx.nj    <- IDX$n[[j]]
    idx.kj    <- IDX$k[[j]]
    idx.skewj <- IDX$skew[[j]]
    idx.ncatj <- IDX$ncat[[j]]

    # `id` is derived from (i, j) so each batch is independent and the ids match
    # the sequential run exactly.
    id    <- (i - 1L) * K + j
    skew  <- as.character(idx.skewj)
    ncat  <- as.character(idx.ncatj)
    n.i   <- n[[idx.nj]]
    k.i   <- k[[idx.kj]]
    model <- models[[idx.modj]]
    seed  <- seeds[[id]]

    print_sep()
    cat(sprintf("i=%i, j=%i, id=%i, total=%i, seed = %i\n", i, j, id, total, seeds[id]))
    print_sep()

    set.seed(seeds[id])

    if (ncat == "Inf") {
      data_i <- sim_cont_data(
        syntax  = model$syntax,
        icc     = model$icc,
        rsd     = model$rsd,
        cluster = "cluster",
        n       = n.i,
        k       = k.i
      )

      ordered <- NULL

    } else {
      thr <- list_thresholds[[skew]][[ncat]]

      data_i <- sim_ord_data(
        syntax  = model$syntax,
        icc     = model$icc,
        rsd     = model$rsd,
        cluster = "cluster",
        n       = n.i,
        k       = k.i,
        thr     = thr
      )

      ordered <- setdiff(colnames(data_i), "cluster")
    }

    print_sep()
    cat(sprintf("Iteration %d/%d:\n", id, total))
    print_sep()

    results.ij <- run_estimators(
      methods  = c("MC-PLS", "Mplus"),
      data_i   = data_i,
      ordered  = ordered,
      model    = model,
      id       = id,
      n.i      = n.i,
      k.i      = k.i,
      skew     = skew,
      ncat     = ncat,
      model.id = idx.modj,
      seed     = seeds[id]
    )

    print(plssem:::plssemParTable(do.call(rbind, unname(results.ij))))
    results.i <- rbind(results.i, do.call(rbind, unname(results.ij)))

    # A failed method has no estimates (see `get_output()`)
    status <- vapply(results.ij, FUN.VALUE = character(1L), FUN = \(r) {
      if (all(is.na(r$est))) "FAILED"
      else if (all(r$admissible)) "admissible"
      else "INADMISSIBLE"
    })

    write_log(sprintf(
      "batch %d, iteration %d/%d (id = %d): model = %s, n = %d, k = %d, ncat = %s, skew = %s | %s",
      i, j, K, id, names(models)[[idx.modj]], n.i, k.i, ncat, skew,
      paste(sprintf("%s: %.1fs, %s", names(results.ij),
                    vapply(results.ij, \(r) r$time[[1L]], numeric(1L)), status),
            collapse = " | ")
    ))
  }

  stamp <- substr(Sys.time(), 1, 16) |>
    stringr::str_replace_all(" ", "-") |>
    stringr::str_replace_all(":", "-")

  filename.sub <- file.path(RESULTS_DIR, sprintf("%s-%s.csv", filePrefix, stamp))
  write.csv(results.i, filename.sub)
  write_log(sprintf("batch %d: finished, wrote %s", i, basename(filename.sub)))

  results.i
}


if (parallel) {
  library(future)
  library(future.apply)

  # Source the project utilities as plain globals so `future` exports them to
  # the workers cleanly, rather than via the devtools::load_all() shadow package
  # (whose name contains a "-" and cannot be attached by name in a worker).
  source(file.path(PROJECT_ROOT, "R", "utils.R"))

  oplan <- plan(multisession, workers = n.workers)
  on.exit(plan(oplan), add = TRUE)

  message(sprintf("Running %d batches across %d workers...", R, n.workers))

  results_list <- future.apply::future_lapply(
    seq_len(R),
    run_batch,
    future.seed     = TRUE,
    future.packages = c("modsem", "plssem", "mvtnorm", "tidyr",
                        "dplyr", "purrr", "stringr", "stats")
  )

} else {
  results_list <- lapply(seq_len(R), run_batch)
}

# Each batch already wrote its own CSV; this is just the in-memory aggregate.
results <- do.call(rbind, results_list)
write_log("Session finished")
