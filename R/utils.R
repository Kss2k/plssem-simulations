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
                       seed     = NULL,
                       true     = NULL,
                       extra    = list()) {
  # `true` can be used to pass the true values directly, as a data.frame with
  # (at least) the columns `par` and `true`. `model` is then passed as is to
  # `func`. Its other columns (e.g., `type`, `level`) are added to the output,
  # and the estimates are matched by `true$key` (if present) or `true$par`.
  # This is needed for multilevel models, which can't be (re)constructed from
  # their parameter table.
  if (is.null(true)) {
    parTable <- modsemify(model)
    parTable2 <- parTable; parTable2$mod <- ""
    syntax.clean <- modsem:::parTableToSyntax(parTable2)
    true <- data.frame(
      par  = paste0(parTable$lhs, parTable$op, parTable$rhs),
      true = as.numeric(parTable$mod)
    )
  } else {
    syntax.clean <- model
  }

  key <- if (is.null(true$key)) true$par else true$key

  output <- data.frame(
    par = true$par,
    est = NA,
    se  = NA,
    true = true$true,
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

  info <- c(true[setdiff(names(true), c("par", "true", "key"))], extra)
  output[names(info)] <- info

  if (!is.null(seed))
    set.seed(seed)

  f.quiet <- purrr::quietly(func)

  time <- system.time({
    tryCatch(
      expr = {
        results <- f.quiet(syntax.clean, data, ...)

        est <- results$result
        idx <- match(key, rownames(est)) # parameters which are missing in `est` are NA
        output$est <- est[idx, "est"]
        output$se  <- est[idx, "se"]
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


is_mlm_syntax <- function(syntax) {
  parTable <- modsemify(syntax)
  any(parTable$lhs == "level" & parTable$op == ":")
}


is_random_slope <- function(mod) startsWith(mod, "rv(")


# Index of the level (0 = before `level: 1`, 1, 2) for each row of a parameter table
level_index <- function(parTable) {
  cumsum(parTable$lhs == "level" & parTable$op == ":")
}


# Get the estimation syntax from a population syntax, by removing the
# (numeric) population values, e.g., `0.4*X` -> `X`. Random slopes are
# specified with duplicated statements, so dulicates are removed as well
get_estimation_syntax <- function(syntax) {
  parTable <- modsemify(syntax)

  key    <- paste(level_index(parTable), parTable$lhs, parTable$op, parTable$rhs)
  is.rv  <- is_random_slope(parTable$mod)
  is.dup <- !is.rv & key %in% key[is.rv]

  parTable <- parTable[!is.dup, , drop = FALSE]
  parTable$mod[!is_random_slope(parTable$mod)] <- ""

  modsem:::parTableToSyntax(parTable)
}


# Split the parameter table of a multilevel syntax into its levels
split_levels <- function(parTable) {
  is.level <- parTable$lhs == "level" & parTable$op == ":"
  idx.l1   <- which(is.level & parTable$rhs %in% c("1", "within"))
  idx.l2   <- which(is.level & parTable$rhs %in% c("2", "between"))

  stopifnot(length(idx.l1) == 1L, length(idx.l2) == 1L, idx.l1 < idx.l2)

  list(
    level.1 = parTable[(idx.l1 + 1L):(idx.l2 - 1L), , drop = FALSE],
    level.2 = parTable[(idx.l2 + 1L):NROW(parTable), , drop = FALSE]
  )
}


# `n` is the sample size within each cluster, and `k` the number of clusters
# (only used for multilevel models). For multilevel models, the parameters at
# each level are standardized (within and between) and `icc` is the share of
# the (total) variance of the indicators at the between level. 
sim_cont_data <- function(syntax, n, k = NULL, cluster = "cluster", icc = NULL, rsd = NULL) {
  parTable <- modsemify(syntax)
  parTable$est <- suppressWarnings(as.numeric(parTable$mod))

  if (!is_mlm_syntax(syntax))
    return(plssem:::simulateDataParTable(parTable, N = n)$ov)

  stopifnot(!is.null(k), !is.null(icc))

  # The structure of the model (e.g., which indicators are at which levels)
  syntax.clean <- get_estimation_syntax(syntax)
  vars <- c(plssem:::getOVs(modsemify(syntax.clean)), cluster)
  dataTemplate <- as.data.frame(
    matrix(0, nrow = 1L, ncol = length(vars), dimnames = list(NULL, vars))
  )

  parsed <- plssem:::parseMultilevelModelArguments(
    syntax  = syntax.clean,
    data    = dataTemplate,
    cluster = cluster
  )

  if (length(icc) == 1L && is.null(names(icc)))
    icc <- stats::setNames(rep(icc, length(parsed$ovs.both)), nm = parsed$ovs.both)

  stopifnot(all(parsed$ovs.both %in% names(icc)))

  # Random slopes. The `rv()` rows are dropped, as the fixed effects are given
  # by the duplicated statements, and the random effects by `rsd`
  rslopes <- parsed$rslopes
  if (NROW(rslopes)) {
    if (is.null(rsd) || !all(rslopes$name %in% names(rsd)))
      stop("`rsd` must be given for all random slopes: ", paste(rslopes$name, collapse = ", "))

    key   <- paste(level_index(parTable), parTable$lhs, parTable$op, parTable$rhs)
    is.rv <- is_random_slope(parTable$mod)
    fixed <- key[!is.rv & !is.na(parTable$est)]

    if (!all(key[is.rv] %in% fixed))
      stop("All random slopes need a (duplicated) fixed effect, e.g., `Y ~ rv(\"s\")*X + 0.4*X`!")

    parTable <- parTable[!is.rv, , drop = FALSE]
    rsd      <- rsd[rslopes$name]

  } else {
    rsd <- NULL
  }

  # The parameters which are not specified (e.g., residual variances)
  # are set such that all variables are standardized at each level
  levels  <- split_levels(parTable)
  asFixed <- \(pt) data.frame(pt[c("lhs", "op", "rhs", "est")], is.free = FALSE)

  parStruct <- list(
    level.1 = asFixed(levels$level.1),
    level.2 = asFixed(levels$level.2),
    icc     = icc[parsed$ovs.both],
    rsd     = rsd
  )

  sim <- plssem:::simulateMultilevelData(
    parStruct    = parStruct,
    parsed       = parsed,
    clusterSizes = rep(n, k),
    times        = 1L,
    seed         = NULL,
    is.hi.ord    = c(FALSE, FALSE),
    full         = c(FALSE, FALSE)
  )

  data <- as.data.frame(sim$ov)[parsed$ovs.all]
  data[[cluster]] <- sim$clusterIdx

  data
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


sim_ord_data <- function(syntax, n, thr = NULL, choose = NULL, cluster = "cluster", ...) {
  cont <- sim_cont_data(syntax = syntax, n = n, cluster = cluster, ...)

  if (is.null(thr)) return(cont)

  # never cut the cluster variable
  if (is.null(choose))
    choose <- setdiff(colnames(cont), cluster)

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
