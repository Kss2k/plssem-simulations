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
                       method = NA,
                       id     = NA,
                       skew   = "",
                       ncat   = "",
                       n      = 0,
                       seed   = NULL) {
  parTable <- modsemify(model)
  parTable2 <- parTable; parTable2$mod <- ""
  syntax.clean <- modsem:::parTableToSyntax(parTable2)

  output <- data.frame(
    par = paste0(parTable$lhs, parTable$op, parTable$rhs),
    est = NA,
    true = as.numeric(parTable$mod),
    method = method,
    id = id,
    seed = seed,
    n = n,
    skew = skew,
    ncat = ncat,
    admissible = FALSE,
    time = NA
  )

  if (!is.null(seed))
    set.seed(seed)

  time <- system.time({
    tryCatch(
      expr = {
        est <- func(syntax.clean, data, ...)
        output$est <- est[output$par]
        output$admissible <- attr(est, "admissible")
      },
      error = \(e) {
        warning(sprintf("%s (%d) failed!, message:\n %s", method, id, e))
        NULL
      }
    )
  })

  output$time <- time

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


print_sep <-  \() cat(strrep("─", options("width")[[1]]), "\n")


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
