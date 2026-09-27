# ──────────────────────────────────────────────────────────────────────────────
# Setup
# ──────────────────────────────────────────────────────────────────────────────
devtools::load_all() # load custom utils

library(mvtnorm)
library(tidyr)
library(dplyr)
library(modsem)        # v.1.0.23
library(plssem)        # v.0.1.5
library(lavaan)        # v.0.7.2  (LSAM)
library(covsim)        # v.1.1.0  (VITA)
library(rvinecopulib)  # sampling from the calibrated vine

.root_info   <- set_project_root()
setwd("mcpls-nnorm")

# Absolute path, so parallel workers (which do not inherit the master's working
# directory) read/write results in the right place regardless of their wd.
PROJECT_ROOT <- .root_info$project.root
RESULTS_DIR <- normalizePath("results", mustWork = FALSE)

# Run settings. Edit these before sourcing/running this script in R.
checkIfExists <- TRUE
R             <- 200L
run.id        <- NULL

parallel  <- TRUE
n.workers <- 8

# ──────────────────────────────────────────────────────────────────────────────
# Model+Parameters
# ──────────────────────────────────────────────────────────────────────────────
# Same population models as mcpls-nlin -- this simulation changes the *data
# generating distribution*, not the model. Only `idx.model` below selects which
# ones are actually run.

# Formula for composite reliability:
cr <- function(...) (sum(c(...))^2)/(sum(c(...))^2 + sum(1-c(...)^2))
cr(0.4, 0.8, 0.8)
#> [1] 0.7194245
cr(0.6, 0.6, 0.8)
#> [1] 0.7092199
cr(0.4, 0.7, 0.9)
#> [1] 0.7220217

models <- c(
  # loadings targeting cr slightly above 0.7
  'X =~ 0.4 * x1 + 0.8 * x2 + 0.8 * x3
   Z =~ 0.6 * z1 + 0.6 * z2 + 0.8 * z3
   Y =~ 0.4 * y1 + 0.7 * y2 + 0.9 * y3
   Y  ~ 0.4 *  X + 0.5 *  Z + 0.3 * X:Z
   X ~~ 0.2 *  Z',
  # reliability 0.9^2, 3 indicators
  'X =~ 0.9 * x1 + 0.9 * x2 + 0.9 * x3
   Z =~ 0.9 * z1 + 0.9 * z2 + 0.9 * z3
   Y =~ 0.9 * y1 + 0.9 * y2 + 0.9 * y3
   Y  ~ 0.4 *  X + 0.5 *  Z + 0.3 * X:Z
   X ~~ 0.2 * Z',
  # reliability 0.8^2, 3 indicators
  'X =~ 0.8 * x1 + 0.8 * x2 + 0.8 * x3
   Z =~ 0.8 * z1 + 0.8 * z2 + 0.8 * z3
   Y =~ 0.8 * y1 + 0.8 * y2 + 0.8 * y3
   Y  ~ 0.4 *  X + 0.5 *  Z + 0.3 * X:Z
   X ~~ 0.2 * Z',
  # reliability 0.5^2, 3 indicators
  'X =~ 0.5 * x1 + 0.5 * x2 + 0.5 * x3
   Z =~ 0.5 * z1 + 0.5 * z2 + 0.5 * z3
   Y =~ 0.5 * y1 + 0.5 * y2 + 0.5 * y3
   Y  ~ 0.4 *  X + 0.5 *  Z + 0.3 * X:Z
   X ~~ 0.2 * Z',
  # reliability 0.9^2, 2 indicators
  'X =~ 0.9 * x1 + 0.9 * x2
   Z =~ 0.9 * z1 + 0.9 * z2
   Y =~ 0.9 * y1 + 0.9 * y2
   Y  ~ 0.4 *  X + 0.5 *  Z + 0.3 * X:Z
   X ~~ 0.2 * Z',
  # reliability 0.8^2, 2 indicators
  'X =~ 0.8 * x1 + 0.8 * x2
   Z =~ 0.8 * z1 + 0.8 * z2
   Y =~ 0.8 * y1 + 0.8 * y2
   Y  ~ 0.4 *  X + 0.5 *  Z + 0.3 * X:Z
   X ~~ 0.2 * Z',
  # reliability 0.5^2, 3 indicators
  'X =~ 0.6 * x1 + 0.5 * x2
   Z =~ 0.6 * z1 + 0.5 * z2
   Y =~ 0.6 * y1 + 0.5 * y2
   Y  ~ 0.4 *  X + 0.5 *  Z + 0.3 * X:Z
   X ~~ 0.2 * Z'
)


# Based on Rhemtulla et al., 2012 and Schubert et al., 2018
# Only `Symmetric` is used here: the distributional conditions below take over
# the role the threshold-asymmetry factor played in mcpls-nlin. The remaining
# entries are kept so `idx.skew` can be widened without further edits.
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


n <- c(300, 1000)

# Set up selection indices which are crossed
idx.model     <- 1
idx.n         <- seq_along(n)
idx.ncat      <- c("2", "3", "5", "7")
idx.skew      <- "Symmetric"
idx.dist.exo  <- NNORM_DIST_EXO  # normal, skewed, uniform
idx.dist.zeta <- NNORM_DIST_ZETA # normal, skewed

# The continuous and the ordinal arm share every factor except `ncat`, which
# only applies once the indicators are discretised. They are therefore built
# separately and stacked, rather than crossed with a dummy `ncat` level.
IDX.continuous <- expand.grid(
  model     = idx.model,
  n         = idx.n,
  ncat      = NA_character_,
  skew      = idx.skew,
  dist.exo  = idx.dist.exo,
  dist.zeta = idx.dist.zeta,
  type      = "continuous",
  stringsAsFactors = FALSE
)

IDX.ordinal <- expand.grid(
  model     = idx.model,
  n         = idx.n,
  ncat      = idx.ncat,
  skew      = idx.skew,
  dist.exo  = idx.dist.exo,
  dist.zeta = idx.dist.zeta,
  type      = "ordinal",
  stringsAsFactors = FALSE
)

IDX <- rbind(IDX.continuous, IDX.ordinal)

# ──────────────────────────────────────────────────────────────────────────────
# Calibration
# ──────────────────────────────────────────────────────────────────────────────
# Two quantities are calibrated once, up front, and then held fixed for every
# replication and every worker:
#
#   vine -- `covsim::vita()` calibrates its pair-copulas by simulation, so the
#           calibrated vine is itself random. Calibrating once keeps the DGP
#           identical across replications instead of drifting between batches.
#
#   psi  -- the structural residual variance, set so that Var(eta) = 1 and the
#           standardized population values stay interpretable. It has to be
#           calibrated *per `dist.exo`*, because Var(X:Z) is a function of the
#           marginal distributions: under right-skewed predictors the product
#           term carries far more variance, so a single psi would leave
#           Var(eta) != 1 in some conditions.
#
# Both are cheap relative to the estimation (a few seconds each) and depend
# only on (model, dist.exo), so there are `length(idx.model) * 3` of them.

VITA_NMAX <- 1e6
VITA_SEED <- 4321L
PSI_N     <- 1e6
PSI_SEED  <- 8765L

calib_key <- function(model.id, dist.exo) paste(model.id, dist.exo, sep = "|")

CALIB <- list()

for (m in idx.model) {
  parTable <- modsemify(models[[m]])
  parTable$est <- as.numeric(parTable$mod)

  info         <- nnorm_struct_info(parTable)
  sigma.target <- nnorm_sigma_target(parTable, xis = info$xis)

  for (dist.exo in idx.dist.exo) {
    message(sprintf("Calibrating vine + psi: model=%d, dist.exo=%s ...", m, dist.exo))

    vine <- calibrate_vine(
      sigma.target = sigma.target,
      dist.exo     = dist.exo,
      Nmax         = VITA_NMAX,
      seed         = VITA_SEED
    )

    psi <- calibrate_psi(
      syntax   = models[[m]],
      vine     = vine,
      dist.exo = dist.exo,
      N        = PSI_N,
      seed     = PSI_SEED
    )

    message(sprintf(
      "  psi = %s",
      paste(sprintf("%s=%.4f", names(psi), psi), collapse = ", ")
    ))

    CALIB[[calib_key(m, dist.exo)]] <- list(vine = vine, psi = psi)
  }
}

# ──────────────────────────────────────────────────────────────────────────────
# Estimators
# ──────────────────────────────────────────────────────────────────────────────

est_pls <- function(model, data, ...) {
  fit <- plssem::pls(model, data, ...)
  par <- plssem::parameter_estimates(fit)

  coef <- cbind(par$est, par$se)
  rownames(coef) <- paste0(par$lhs, par$op, par$rhs)
  colnames(coef) <- c("est", "se")

  par.admissible <- checkIfParTableIsAdmissible(par)
  fit.admissible <- is_admissible(fit) # fit.admissible should be sufficient
                                       # but we check both (just in case)

  if (!par.admissible && fit.admissible) {
    warning("pars are inadmissible! But fit says it's admissible!")
  } else if (!par.admissible && !fit.admissible) {
    warning("both pars and fit are inadmissible!")
  } else if (par.admissible && !fit.admissible) {
    warning("pars are admissible! But fit says it's inadmissible!")
  }

  attr(coef, "admissible") <- fit.admissible && par.admissible

  coef
}


est_mplus <- function(model, data, ...) {
  fit <- modsem::modsem_mplus(model, data, ...)

  mod <- modsemify(model)
  par <- modsem::standardized_estimates(fit)

  # Mplus is case insenstitive, so we have to account for that
  vars0 <- union(mod$lhs, mod$rhs)
  vars1 <- union(par$lhs, par$rhs)

  # Create mapping from upper to lower
  is.upper <- tolower(vars1) %in% vars0
  mapping <- stats::setNames(vars1, nm = vars1)
  mapping[is.upper] <- tolower(mapping[is.upper])

  # Map upper to lower
  par$lhs <- mapping[par$lhs]
  par$rhs <- mapping[par$rhs]

  coef <- cbind(par$est, par$std.error)
  rownames(coef) <- paste0(par$lhs, par$op, par$rhs)
  colnames(coef) <- c("est", "se")

  attr(coef, "admissible") <- checkIfParTableIsAdmissible(par)

  coef
}


# Local structural-after-measurement, via lavaan's `sam()` with its defaults
# (`sam_method = "local"`, `se = "twostep"`, `mm_args = list(bounds =
# "wide.zerovar")`).
#
# `sam()` cannot combine latent interactions with the correlation structures it
# uses for ordered indicators ("SAM + lv interactions do not work (yet) if
# correlation structures are used"), so ordinal indicators are handed to it as
# continuous -- the same treatment PLS and PLSc give them.
#
# The population values are standardized, so the estimates are read off
# `standardizedSolution()`, as for Mplus.
est_lsam <- function(model, data, ...) {
  fit <- lavaan::sam(model, data = data, ...)
  par <- lavaan::standardizedSolution(fit)

  coef <- cbind(par$est.std, par$se)
  rownames(coef) <- paste0(par$lhs, par$op, par$rhs)
  colnames(coef) <- c("est", "se")

  par$est <- par$est.std
  attr(coef, "admissible") <- checkIfParTableIsAdmissible(par)

  coef
}


# Column order of the results CSV. The first block matches mcpls-nlin; `type`,
# `dist.exo` and `dist.zeta` are the new design factors.
NNORM_COLS <- c(
  "par", "est", "se", "true", "method", "id", "seed", "n", "type",
  "dist.exo", "dist.zeta", "skew", "ncat", "model.id", "admissible", "time"
)


# `get_output()` is shared with mcpls-nlin, so rather than widening its
# signature the new design factors are attached here.
get_output_nnorm <- function(..., type, dist.exo, dist.zeta) {
  out <- get_output(...)

  out$type      <- type
  out$dist.exo  <- dist.exo
  out$dist.zeta <- dist.zeta

  out[NNORM_COLS]
}


run_estimators <- function(methods,
                           data_i,
                           ordered,
                           model,
                           id,
                           n.i,
                           type,
                           dist.exo,
                           dist.zeta,
                           skew,
                           ncat,
                           model.id,
                           seed) {
  methods <- match.arg(
    methods,
    choices = c("MC-PLSc", "PLSc", "PLS", "Mplus", "LSAM"),
    several.ok = TRUE
  )

  is.ordinal <- type == "ordinal"

  # Metadata carried on every row, identical across methods.
  meta <- list(
    data      = data_i,
    model     = model,
    id        = id,
    n         = n.i,
    type      = type,
    dist.exo  = dist.exo,
    dist.zeta = dist.zeta,
    skew      = skew,
    ncat      = ncat,
    model.id  = model.id,
    seed      = seed
  )

  run_one <- function(method) {
    args <- switch(
      method,

      # Monte-Carlo consistent PLS. With ordered indicators `pls()` selects the
      # MC algorithm itself; on continuous data it has to be asked for.
      "MC-PLSc" = if (is.ordinal) {
        list(func = est_pls, ordered = ordered, bootstrap = TRUE, boot.R = 500)
      } else {
        list(func = est_pls, mcpls = TRUE, bootstrap = TRUE, boot.R = 500)
      },

      "PLSc" = list(
        func       = est_pls,
        bootstrap  = TRUE,
        boot.R     = 500,
        consistent = TRUE
      ),

      "PLS" = list(
        func       = est_pls,
        bootstrap  = TRUE,
        boot.R     = 500,
        consistent = FALSE
      ),

      "Mplus" = list(
        func        = est_mplus,
        cleanup     = TRUE,
        estimator   = "mlr",
        processors  = 1,
        categorical = if (is.ordinal) ordered else NULL
      ),

      # lavaan's sam() defaults; ordinal indicators treated as continuous
      "LSAM" = list(func = est_lsam)
    )

    do.call(get_output_nnorm, c(args, meta, list(method = method)))
  }

  stats::setNames(lapply(methods, run_one), methods)
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
# v2-test is for testing. The other run.ids (see below) specify
# different computers running parallel simulations. They of course
# need different seeds, to generate unique results. These differ from the
# mcpls-nlin seeds, so the two studies do not share draws.

DEFAULT_SEED_IDX <- 5
LOCAL_SEEDS <- c(
  "v2-test"     = 7715284,
  "v2-vivo"     = 3092648,
  "v2-tuf"      = 6481073,
  "v2-promax"   = 1837520,
  "v2-lovelace" = 2139260
)

if (is.null(run.id)) {
  cat("What run.id do you want to use? Available:\n")
  print(names(LOCAL_SEEDS))
  run.id.idx <- as.integer(readLines(n=1))

  if (!length(run.id.idx) || is.na(run.id.idx)) run.id.idx <- DEFAULT_SEED_IDX
  run.id <- names(LOCAL_SEEDS)[[run.id.idx]]
}

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

    results.i <- read.csv(file.path(RESULTS_DIR, files[which(match)[[1L]]]))
    results.i <- results.i[-1] # drop rownames
    return(results.i)
  }

  for (j in seq_len(NROW(IDX))) {

    idx.modj   <- IDX$model[[j]]
    idx.nj     <- IDX$n[[j]]
    typej      <- IDX$type[[j]]
    skewj      <- IDX$skew[[j]]
    ncatj      <- IDX$ncat[[j]]
    dist.exoj  <- IDX$dist.exo[[j]]
    dist.zetaj <- IDX$dist.zeta[[j]]

    # `id` is derived from (i, j) so each batch is independent and the ids match
    # the sequential run exactly.
    id    <- (i - 1L) * K + j
    n.i   <- n[[idx.nj]]
    model <- models[[idx.modj]]
    seed  <- seeds[[id]]
    calib <- CALIB[[calib_key(idx.modj, dist.exoj)]]

    print_sep()
    cat(sprintf(
      "i=%i, j=%i, id=%i, total=%i, seed=%i, type=%s, dist.exo=%s, dist.zeta=%s\n",
      i, j, id, total, seeds[id], typej, dist.exoj, dist.zetaj
    ))
    print_sep()

    set.seed(seeds[id])

    # `thr = NULL` leaves the indicators continuous
    thr <- if (typej == "ordinal") list_thresholds[[skewj]][[ncatj]] else NULL

    data_i <- sim_nnorm_data(
      syntax    = model,
      n         = n.i,
      vine      = calib$vine,
      psi       = calib$psi,
      dist.exo  = dist.exoj,
      dist.zeta = dist.zetaj,
      thr       = thr
    )

    ordered <- if (typej == "ordinal") colnames(data_i) else NULL

    print_sep()
    cat(sprintf("Iteration %d/%d:\n", id, total))
    print_sep()

    results.ij <- run_estimators(
      methods   = c("MC-PLSc", "PLSc", "PLS", "Mplus", "LSAM"),
      data_i    = data_i,
      ordered   = ordered,
      model     = model,
      id        = id,
      n.i       = n.i,
      type      = typej,
      dist.exo  = dist.exoj,
      dist.zeta = dist.zetaj,
      skew      = skewj,
      ncat      = ncatj,
      model.id  = idx.modj,
      seed      = seeds[id]
    )

    print(plssem:::plssemParTable(do.call(rbind, unname(results.ij))))
    results.i <- rbind(results.i, do.call(rbind, unname(results.ij)))
  }

  stamp <- substr(Sys.time(), 1, 16) |>
    stringr::str_replace_all(" ", "-") |>
    stringr::str_replace_all(":", "-")

  filename.sub <- file.path(RESULTS_DIR, sprintf("%s-%s.csv", filePrefix, stamp))
  write.csv(results.i, filename.sub)

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

  message(sprintf(
    "%s: Running %d batches across %d workers...",
    run.id, R, n.workers
  ))

  results_list <- future.apply::future_lapply(
    seq_len(R),
    run_batch,
    future.seed     = TRUE,
    future.packages = c("modsem", "plssem", "lavaan", "covsim", "rvinecopulib",
                        "mvtnorm", "tidyr", "dplyr", "purrr", "stringr", "stats")
  )

} else {
  results_list <- lapply(seq_len(R), run_batch)
}

# Each batch already wrote its own CSV; this is just the in-memory aggregate.
results <- do.call(rbind, results_list)
