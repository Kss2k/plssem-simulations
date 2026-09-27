set_project_root <- function(target = "DESCRIPTION", max.iter=100) {
  user_root    <- normalizePath("~")
  original_dir <- getwd()

  found <- FALSE
  i     <- 0

  while (!found && i < max.iter) {
    found <- target %in% list.files()
  
    if ((i <- i + 1) > max.iter || found || 
      getwd() == user_root || getwd() == "/") break
    
    setwd("..")
  }

  if (!found) {
    setwd(original_dir)
    warning("Project root not found!")
    project_root <- NA 

  } else {
    message("Project root found!")
    project_root <- getwd() 
  } 

  invisible(list(
    success = found,
    project.root = project_root,
    original.dir = original_dir,
    user.root = user_root,
    target = target,
    iter = i
  ))
}


get_output <- function(func,
                       model,
                       data,
                       ...,
                       method   = NA,
                       id       = NA,
                       skew     = "",
                       ncat     = "",
                       n        = 0,
                       model.id = 0,
                       seed     = NULL) {
  parTable <- modsemify(model)
  parTable2 <- parTable; parTable2$mod <- ""
  syntax.clean <- modsem:::parTableToSyntax(parTable2)

  output <- data.frame(
    par = paste0(parTable$lhs, parTable$op, parTable$rhs),
    est = NA,
    se  = NA,
    true = as.numeric(parTable$mod),
    method = method,
    id = id,
    seed = seed,
    n = n,
    skew = skew,
    ncat = ncat,
    model.id = model.id,
    admissible = FALSE,
    time = NA
  )

  if (!is.null(seed))
    set.seed(seed)

  f.quiet <- purrr::quietly(func)

  time <- system.time({
    tryCatch(
      expr = {
        results <- f.quiet(syntax.clean, data, ...)

        est <- results$result
        output$est <- est[output$par, "est"]
        output$se  <- est[output$par, "se"]
        output$admissible <- attr(est, "admissible")

        if (length(results$warnings)) {
          msg <- paste0(results$warnings, collapse = ";")
          plssem:::pls_msg_warn_immediate(
            sprintf("method=%s, id=%i, message(s)=%s", method, id, msg)
          )
        }
      },
      error = \(e) {
        warning(sprintf("%s (%d) failed!, message:\n %s", method, id, e))
        NULL
      }
    )
  })

  output$time <- time[["elapsed"]]

  class(output) <- c("simoutput", "data.frame")
  output
}


print.simoutput <- function(x, ...) {
  cat(sprintf(
    "ID: %i, Method: %s, Elapsed: %s Skew: %s, NCAT: %d, N: %d\n",
    unique(x$id), unique(x$method), capture.output(x$elapsed),
    unique(x$skew), unique(x$ncat), unique(x$n)
  ))

  print(as.data.frame(x))
}


print_sep <- \() cat(strrep("─", options("width")[[1]]), "\n")


sim_cont_data <- function(syntax, n = n) {
  
  parTable <- modsemify(syntax)
  parTable$est <- as.numeric(parTable$mod)

  plssem:::simulateDataParTable(parTable, N = n)$ov
}


cut_data <- function(data, thr, choose = NULL) {
  standardize <- \(x) (x - mean(x)) / sd(x)

  if (is.null(choose))
    choose <- colnames(data)

  for (i in seq_along(choose)) {
    var <- choose[[i]]
    x <- data[[var]]
    breaks <- c(-Inf, thr, Inf)
    y <- cut(standardize(x), breaks = breaks)
    y <- as.integer(as.ordered(as.integer(y)))

    data[[var]] <- y
    z <- rep(NA_real_, length(y))

  }

  data
}


sim_ord_data <- function(syntax, n, thr = NULL, choose = NULL) {
  cont <- sim_cont_data(syntax = syntax, n = n)

  if (is.null(thr)) return(cont)

  cut_data(cont, thr = thr, choose = choose)
}


checkIfParTableIsAdmissible <- function(parTable, n = 5000) {
  # use the simulateDataParTable utility function to check if the
  # parameter estimates are able to simulate an admissible parTable
  tryCatch({
    sim <- plssem:::simulateDataParTable(parTable, N = n)
    sim$is.admissible
  }, error = function(e) {
    warning("Failed to simulate data to check admissibility!\n",
            "Message: ", conditionMessage(e))
    FALSE
  })
}


# ──────────────────────────────────────────────────────────────────────────────
# Non-normal data generation (covsim/VITA)
# ──────────────────────────────────────────────────────────────────────────────
# Follows the non-normality specification of Vieira, Slupphaug & Rosseel: the
# exogenous predictors are drawn through a calibrated vine copula (VITA), and
# the structural disturbances are drawn independently, with the two forming a
# factorial design. Measurement errors stay normal -- the specification covers
# the structural part of the model only.

NNORM_DIST_EXO  <- c("normal", "skewed", "uniform")
NNORM_DIST_ZETA <- c("normal", "skewed")


# Marginal specifications passed to `covsim::vita()`. Every condition has unit
# variance, so the margins line up with the standardized latent variables of
# the population models without any rescaling:
#   normal  : N(0, 1)                 skew 0, excess kurtosis  0
#   skewed  : gamma(1, 1)             skew 2, excess kurtosis  6
#   uniform : U(-sqrt(3), sqrt(3))    skew 0, excess kurtosis -1.2
nnorm_margins <- function(dist.exo, p) {
  margin <- switch(
    dist.exo,
    normal  = list(distr = "norm",  mean = 0, sd = 1),
    skewed  = list(distr = "gamma", shape = 1, rate = 1),
    uniform = list(distr = "unif",  min = -sqrt(3), max = sqrt(3)),
    stop("Unknown dist.exo: ", dist.exo, call. = FALSE)
  )

  rep(list(margin), p)
}


# Population (not sample) mean of each marginal. Subtracted after sampling so
# the predictors are centred at zero in the population, keeping the conditions
# comparable while leaving the sample means free to fluctuate.
nnorm_mean <- function(dist.exo) {
  switch(
    dist.exo,
    normal  = 0,
    skewed  = 1, # gamma(shape, rate) has mean shape/rate
    uniform = 0,
    stop("Unknown dist.exo: ", dist.exo, call. = FALSE)
  )
}


# `vita()` calibrates its pair-copulas by simulation, so the calibrated vine is
# itself random. We calibrate once per condition under a fixed seed and reuse
# the vine for every replication, which keeps the DGP identical across
# replications and across parallel workers.
#
# `family_set = "gauss"` in every condition, so the distribution factor varies
# the marginals *only* and the `normal` condition is exactly multivariate
# normal (reproducing the DGP used in mcpls-nlin). Note that `vita()` defaults
# to trying "clayton" first, which would otherwise leave even the `normal`
# condition non-normal in its dependence structure.
#
# `cores = 1` matters: the default is `parallel::detectCores()`, which would
# oversubscribe the machine when called from inside a `future` worker.
calibrate_vine <- function(sigma.target, dist.exo, Nmax = 1e6, seed = 4321L) {
  margins <- nnorm_margins(dist.exo, p = NCOL(sigma.target))

  set.seed(seed)
  covsim::vita(
    margins,
    sigma.target = sigma.target,
    family_set   = "gauss",
    Nmax         = Nmax,
    cores        = 1L,
    verbose      = FALSE
  )
}


# Target covariance matrix of the exogenous latent variables, read off the
# `~~` rows of the parameter table (unit variances, as the models are
# standardized).
nnorm_sigma_target <- function(parTable, xis) {
  mat <- diag(1, length(xis))
  dimnames(mat) <- list(xis, xis)

  for (i in seq_along(xis)) for (j in seq_len(i - 1L)) {
    cond <- (
      (parTable$op == "~~" & parTable$lhs == xis[[i]] & parTable$rhs == xis[[j]]) |
      (parTable$op == "~~" & parTable$lhs == xis[[j]] & parTable$rhs == xis[[i]])
    )

    cov.ij <- if (any(cond)) parTable[cond, "est"][[1L]] else 0
    mat[i, j] <- mat[j, i] <- cov.ij
  }

  mat
}


draw_exo <- function(vine, dist.exo, n, xis) {
  x <- rvinecopulib::rvine(n = n, vine = vine)
  x <- x - nnorm_mean(dist.exo)

  colnames(x) <- xis
  as.data.frame(x)
}


# Structural disturbance. Both options have mean 0 and variance `psi`:
#   normal : N(0, psi)
#   skewed : a centred, scaled exponential -- rate 1/sqrt(psi), shifted by
#            -sqrt(psi) -- with skew 2 and excess kurtosis 6, i.e. the same
#            base distribution as the right-skewed exogenous predictors.
draw_zeta <- function(n, psi, dist.zeta) {
  switch(
    dist.zeta,
    normal = stats::rnorm(n, mean = 0, sd = sqrt(psi)),
    skewed = stats::rexp(n, rate = 1 / sqrt(psi)) - sqrt(psi),
    stop("Unknown dist.zeta: ", dist.zeta, call. = FALSE)
  )
}


nnorm_struct_info <- function(parTable) {
  struct <- parTable[parTable$op == "~", , drop = FALSE]
  etas   <- unique(struct$lhs)
  base   <- unique(unlist(strsplit(unique(struct$rhs), ":", fixed = TRUE)))

  list(struct = struct, etas = etas, xis = setdiff(base, etas))
}


# Build the endogenous latent variables from the exogenous draws. When `psi` is
# NULL the residual variances are *calibrated* instead of used: each is set to
# 1 - Var(deterministic part), so that the latent variable has unit variance
# and the standardized population values stay interpretable. Product terms are
# mean-centred, as in `plssem:::simulateDataParTable()`.
build_latents <- function(Xi, info, psi = NULL, dist.zeta = "normal") {
  calibrate <- is.null(psi)
  n <- NROW(Xi)

  if (calibrate) {
    psi       <- stats::setNames(rep(NA_real_, length(info$etas)), info$etas)
    dist.zeta <- "normal" # the shape of zeta does not affect its variance
  }

  for (eta in info$etas) {
    rows <- info$struct[info$struct$lhs == eta, , drop = FALSE]

    for (term in rows$rhs[grepl(":", rows$rhs, fixed = TRUE)]) {
      if (!is.null(Xi[[term]])) next

      elems      <- strsplit(term, ":", fixed = TRUE)[[1L]]
      Xi[[term]] <- Reduce(`*`, Xi[elems])
      Xi[[term]] <- Xi[[term]] - mean(Xi[[term]])
    }

    vals <- as.vector(as.matrix(Xi[rows$rhs]) %*% rows$est)

    if (calibrate) {
      psi[[eta]] <- 1 - stats::var(vals)

      if (psi[[eta]] <= 0)
        warning(sprintf(
          "Calibrated residual variance for %s is not positive (%.4f)!",
          eta, psi[[eta]]
        ))
    }

    Xi[[eta]] <- vals + draw_zeta(n, psi = max(psi[[eta]], 0), dist.zeta = dist.zeta)
  }

  list(Xi = Xi, psi = psi)
}


# Residual variances, calibrated once per (model, dist.exo) from a large
# population draw and then held fixed across every replication -- the fixed,
# calibrated `psi` of the specification, rather than a quantity re-derived from
# each sample. Depends on `dist.exo` because Var(X:Z) is a function of the
# marginals, so a single `psi` would leave Var(eta) != 1 in some conditions.
calibrate_psi <- function(syntax, vine, dist.exo, N = 1e6, seed = 8765L) {
  parTable <- modsemify(syntax)
  parTable$est <- as.numeric(parTable$mod)

  info <- nnorm_struct_info(parTable)

  set.seed(seed)
  Xi <- draw_exo(vine, dist.exo, n = N, xis = info$xis)

  build_latents(Xi, info = info, psi = NULL)$psi
}


# Generate one data set. `thr = NULL` returns continuous indicators; otherwise
# the indicators are discretised with `cut_data()`, exactly as in mcpls-nlin.
sim_nnorm_data <- function(syntax,
                           n,
                           vine,
                           psi,
                           dist.exo,
                           dist.zeta,
                           thr    = NULL,
                           choose = NULL) {
  parTable <- modsemify(syntax)
  parTable$est <- as.numeric(parTable$mod)

  info <- nnorm_struct_info(parTable)

  Xi <- draw_exo(vine, dist.exo, n = n, xis = info$xis)
  Xi <- build_latents(Xi, info = info, psi = psi, dist.zeta = dist.zeta)$Xi

  # Measurement model. Errors are normal in every condition, and scaled so the
  # indicators have unit variance given standardized latent variables.
  measr <- parTable[parTable$op == "=~", , drop = FALSE]
  Inds  <- list()

  for (i in seq_len(NROW(measr))) {
    lv     <- measr$lhs[[i]]
    ind    <- measr$rhs[[i]]
    lambda <- measr$est[[i]]

    Inds[[ind]] <- lambda * Xi[[lv]] +
      stats::rnorm(n, mean = 0, sd = sqrt(max(1 - lambda^2, 0)))
  }

  data <- as.data.frame(Inds)

  if (is.null(thr)) return(data)

  cut_data(data, thr = thr, choose = choose)
}
