# ──────────────────────────────────────────────────────────────────────────────
# Setup
# ──────────────────────────────────────────────────────────────────────────────
devtools::load_all() # load custom utils

library(mvtnorm)
library(tidyr)
library(dplyr)
library(modsem) # v.1.0.21
library(plssem) # v.0.1.3

.root_info   <- set_project_root()
setwd("mcpls-nlin")

# Absolute path, so parallel workers (which do not inherit the master's working
# directory) read/write results in the right place regardless of their wd.
PROJECT_ROOT <- .root_info$project.root
RESULTS_DIR <- normalizePath("results", mustWork = FALSE)

# ──────────────────────────────────────────────────────────────────────────────
# Model+Parameters
# ──────────────────────────────────────────────────────────────────────────────

models <- c(
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


n <- c(200, 500, 1000)

# Set up selection indices which are crossed
# idx.model <- c(5, 6) # two indicators, rel = 0.8^2 and 0.6^2
idx.model <- c(2, 3) # two indicators, rel = 0.8^2 and 0.6^2
idx.n     <- seq_along(n)
idx.ncat  <- c("2", "3", "5", "7")
idx.skew  <- c("Symmetric", "Moderate", "Extreme", "Alt.Mod", "Alt.Ext")

IDX <- expand.grid(
  model = idx.model,
  n     = idx.n,
  ncat  = idx.ncat,
  skew  = idx.skew
)

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
                                       # but we check both (just in case) # fit.admissible should be sufficient

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


# ──────────────────────────────────────────────────────────────────────────────
# Run Simulation
# ──────────────────────────────────────────────────────────────────────────────
checkIfExists <- TRUE
R             <- 200L
K             <- NROW(IDX)
total         <- R * K
run.id        <- NULL


# ── Parallelisation ───────────────────────────────────────────────────────────
# Set `parallel <- TRUE` to evaluate the R batches concurrently with the
# `future` package. Each batch (one value of the outer `i` index) is fully
# independent: it derives its iteration ids from `i`, sets its own per-iteration
# seeds, and writes its own results CSV. The batches can therefore run in any
# order / in parallel and still reproduce the sequential output exactly.
parallel  <- TRUE
n.workers <- 4

# The run.id specifies the circumstance the script is running under
# v0-test is for testing. The other run.ids (see below) specify
# different computers running parallel simulations. They of course
# need different seeds, to generate unique results.

LOCAL_SEEDS <- c(
  "v0-test" = 5340956,
  "v0-vivo" = 9144416,
  "v0-tuf"  = 1210967
)

if (is.null(run.id)) {
  cat("What run.id do you want to use? Available:\n")
  print(names(LOCAL_SEEDS))
  run.id.idx <- as.integer(readLines(n=1))
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
  match <- startsWith(files, filePrefix)

  if (checkIfExists && any(match)) {
    message(sprintf(
      "Skipping iteration batch %d, as it has already been run...", i)
    )

    results.i <- read.csv(file.path(RESULTS_DIR, files[which(match)[[1L]]]))
    results.i <- results.i[-1] # drop rownames
    return(results.i)
  }

  for (j in seq_len(NROW(IDX))) {

    idx.modj  <- IDX$model[[j]]
    idx.nj    <- IDX$n[[j]]
    idx.skewj <- IDX$skew[[j]]
    idx.ncatj <- IDX$ncat[[j]]

    # `id` is derived from (i, j) so each batch is independent and the ids match
    # the sequential run exactly.
    id    <- (i - 1L) * K + j
    skew  <- names(list_thresholds)[[idx.skewj]]
    ncat  <- names(list_thresholds[[idx.skewj]][idx.ncatj]) # could just write ncat <- idx.ncatj
    n.i   <- n[[idx.nj]]
    model <- models[[idx.modj]]
    seed  <- seeds[[id]]

    print_sep()
    cat(sprintf("i=%i, j=%i, id=%i, total=%i, seed = %i\n", i, j, id, total, seeds[id]))
    print_sep()

    set.seed(seeds[id])
    thr <- list_thresholds[[skew]][[ncat]]
    data_i <- sim_ord_data(syntax = model, thr = thr, n = n.i)
    ordered <- colnames(data_i)

    print_sep()
    cat(sprintf("Iteration %d/%d:\n", id, total))
    print_sep()

    results.ij <- list(
      mcpls = get_output(
        func      = est_pls,
        data      = data_i,
        ordered   = ordered,
        bootstrap = TRUE,
        boot.R    = 500,
        model     = model,
        method    = "MC-OrdPLSc",
        id        = id,
        n         = n.i,
        skew      = skew,
        ncat      = ncat,
        seed      = seeds[id]
      ),

      plsc = get_output(
        func       = est_pls,
        data       = data_i,
        bootstrap  = TRUE,
        boot.R     = 500,
        model      = model,
        method     = "PLSc",
        consistent = TRUE,
        n          = n.i,
        id         = id,
        skew       = skew,
        ncat       = ncat,
        seed       = seeds[id]
      ),

      pls = get_output(
        func       = est_pls,
        data       = data_i,
        model      = model,
        bootstrap  = TRUE,
        boot.R     = 500,
        method     = "PLS",
        consistent = FALSE,
        n          = n.i,
        id         = id,
        skew       = skew,
        ncat       = ncat,
        seed       = seeds[id]
      ),

      mplus = get_output(
        func        = est_mplus,
        data        = data_i,
        model       = model,
        method      = "Mplus",
        processors  = 8,
        categorical = ordered,
        n           = n.i,
        id          = id,
        skew        = skew,
        ncat        = ncat,
        seed        = seeds[id]
      )
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
