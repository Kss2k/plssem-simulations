devtools::load_all()
set_project_root()

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(scales)

testfiles <- FALSE
rdir <- "mcpls-nnorm/results/"
files <- dir(rdir)
files <- files[endsWith(files, ".csv")]
if (!testfiles) files <- files[!grepl("test", files) & !grepl("extra", files)]

paths <- paste0(rdir, files)

# results-<run.id>-<batch>-<timestamp>.csv
ids <- stats::setNames(paste(
  stringr::str_extract(files, "^results-(v2-[a-z]+)-[0-9]+-", group = 1),
  stringr::str_extract(files, "^results-v2-[a-z]+-([0-9]+)-", group = 1),
  sep = "-"
), nm = files)

read <- function(path) {
  name <- last(stringr::str_split_1(path, "/"))
  df <- read.csv(path)
  id <- ids[[name]]
  df$id <- paste0(id, "-", df$id)
  df
}


methods_drop <- "LSAM" # for now
methods_ordered <- c("PLS", "PLSc", "MC-PLSc", "Mplus")

# The two new design factors. `dist.exo` is the distribution of the exogenous
# predictors (drawn through the calibrated vine) and `dist.zeta` that of the
# structural disturbance. They are crossed with the threshold-asymmetry factor
# `skew` carried over from mcpls-nlin: the ordinal arm varies both, since
# symmetric thresholds alone would artificially favour PLSc and LSAM. The
# continuous arm has no thresholds, so `skew` is constant ("Symmetric") there.
dist_exo_ordered  <- c("normal", "skewed", "uniform")
dist_zeta_ordered <- c("normal", "skewed")
skew_ordered      <- c("Symmetric", "Moderate", "Extreme", "Alt.Mod", "Alt.Ext")

df <- do.call(rbind, lapply(paths, read)) |>
  filter(!method %in% methods_drop) |>
  mutate(
    bias = est - true,
    method = factor(method,
      levels = methods_ordered,
      labels = methods_ordered
    ),
    dist.exo  = factor(dist.exo,  levels = dist_exo_ordered),
    dist.zeta = factor(dist.zeta, levels = dist_zeta_ordered),
    # `intersect` so this keeps working if `idx.skew` is widened again
    skew      = factor(skew, levels = intersect(skew_ordered, unique(skew))),
  ) |>
  group_by(id, method) |>
  mutate(
    parcombo = paste0(paste0(par, "=", true), collapse = ","),
    admissible.se = all(admissible) & !any(is.na(se) | se > 1) # check SEs when checking admissiblity
  )

# We try to split the simulations into data type, sample size, and model
# parameter combos. Within each we look at the performance.
simsplit <- expand.grid(
  type = sort(unique(df$type)),
  n = sort(unique(df$n)),
  model.id = sort(unique(df$model.id)),
  drop.inadmissible = c(TRUE, FALSE),
  stringsAsFactors = FALSE
)

par2tex <- list(
  `Y~X`   = "gamma[1]",
  `Y~Z`   = "gamma[2]",
  `Y~X:Z` = "gamma[3]",
  `Y=~y1` = "lambda[7]",
  `Y=~y2` = "lambda[8]",
  `Y=~y3` = "lambda[9]"
)

# The continuous arm has no `ncat`, so the distribution of the predictors takes
# over the x-axis there and only `dist.zeta` is facetted. In the ordinal arm the
# layout matches mcpls-nlin (categories on x) with the 3 x 2 distributional
# factorial down the rows, and the threshold conditions sitting next to the
# parameter across the columns -- so Symmetric and Extreme end up adjacent and
# can be read off against each other directly.
x_var  <- function(type) if (type == "ordinal") "ncat" else "dist.exo"
x_lab  <- function(type) if (type == "ordinal") "Categories" else "Predictor distribution"
row_facets <- function(type) {
  if (type == "ordinal") vars(dist.exo, dist.zeta) else vars(dist.zeta)
}
# Column facets for the plots that have a parameter dimension, and for those
# that do not (inadmissibility, computation time).
col_facets <- function(type) {
  if (type == "ordinal") vars(par, skew) else vars(par)
}
plain_col_facets <- function(type) {
  if (type == "ordinal") vars(skew) else NULL
}
# `par` holds plotmath ("gamma[1]"); `skew` holds ordinary labels.
facet_labeller <- labeller(par = label_parsed, .default = label_value)

# Computation times are only comparable within one machine, so the timing plots
# are restricted to a single run.id -- whichever contributed most rows, rather
# than a hardcoded name that goes stale when the run moves to another machine.
# `id` is "<run.id>-<batch>-<iteration>", so drop the two trailing numbers.
df$run.id <- sub("-[0-9]+-[0-9]+$", "", df$id)
timing.runid <- names(sort(table(df$run.id), decreasing = TRUE))[[1L]]
message("Timing plots restricted to run.id: ", timing.runid,
        " (of: ", paste(sort(unique(df$run.id)), collapse = ", "), ")")

# Count inadmissibles
admissible <- group_by(df, id, method, model.id, type, ncat, skew, dist.exo, dist.zeta, n) |>
  summarize(admissible = unique(admissible)) |>
  group_by(method, model.id, type, ncat, skew, dist.exo, dist.zeta, n) |>
  summarize(nruns = length(admissible),
            ninadmissible = sum(!admissible),
            pinadmissible = sum(!admissible)/length(admissible))

print(admissible, n = 500)


EMPTY_LIST <- vector("list", NROW(simsplit))

plots_nnorm_inadmissible      <- EMPTY_LIST
plots_nnorm_time              <- EMPTY_LIST
plots_nnorm_bias_l1_l2        <- EMPTY_LIST
plots_nnorm_bias_b1           <- EMPTY_LIST
plots_nnorm_bias_b1_b2        <- EMPTY_LIST
plots_nnorm_bias_b2           <- EMPTY_LIST
plots_nnorm_bias_b3           <- EMPTY_LIST
plots_nnorm_se_sd_ratio_b1    <- EMPTY_LIST
plots_nnorm_se_sd_ratio_b2    <- EMPTY_LIST
plots_nnorm_se_sd_ratio_b3    <- EMPTY_LIST
plots_nnorm_se_sd_ratio_b1_b2 <- EMPTY_LIST
plots_nnorm_se_sd_b1          <- EMPTY_LIST
plots_nnorm_se_sd_b2          <- EMPTY_LIST
plots_nnorm_se_sd_b3          <- EMPTY_LIST

for (i in seq_len(NROW(simsplit))) suppressMessages({
  cat(sprintf("%i...\n", i))
  # ----------------------------------------------------------------------------
  # Simulation settings
  # ----------------------------------------------------------------------------

  type.i <- simsplit$type[[i]]
  n.i <- simsplit$n[[i]]
  model.i <- simsplit$model.id[[i]]
  drop.inadmissible <- simsplit$drop.inadmissible[[i]]

  xv <- x_var(type.i)
  xl <- x_lab(type.i)
  rf <- row_facets(type.i)
  cf <- col_facets(type.i)
  pf <- plain_col_facets(type.i)

  # ----------------------------------------------------------------------------
  # Inadmissible Solutions
  # ----------------------------------------------------------------------------

  dodge <- position_dodge(width = 0.9)
  pinadmissible <- admissible |>
    filter(n == n.i, model.id == model.i, type == type.i) |>
    mutate(xvar = as.factor(.data[[xv]])) |>
    ggplot(aes(x = xvar, y = pinadmissible, colour = method, fill = method)) +
    geom_col(alpha = 0.2, position = dodge) +
    facet_grid(rows = rf, cols = pf, scales = "fixed") +
    # coord_cartesian(ylim = c(0, 1)) +
    scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
    # ggtitle(sprintf("Percentage inadmissible solutions (n=%i) model %d", n.i, model.i)) +
    ylab("Percentage inadmissible solutions") +
    xlab(xl) +
    theme_bw()

  pinadmissible2 <-
    admissible |> mutate(
      xvar = as.factor(.data[[xv]]),
      pinadmissible = 100 * pinadmissible,
      pinadmissible.scaled = (pinadmissible - max(pinadmissible))^(1/8)
    ) |>
    filter(n == n.i, model.id == model.i, type == type.i) |>
    ggplot(aes(x = xvar, y = dist.zeta)) +
    geom_tile(aes(fill=pinadmissible)) +
    geom_text(aes(label=paste0(round(pinadmissible,1), "%"))) +
    facet_grid(rows = vars(method), cols = pf) +
    ylab("Disturbance distribution") +
    xlab(xl) +
    scale_fill_gradient(low = "white", high = "red") +
    theme_bw()


  if (drop.inadmissible) {
    ids.is.admissible <- group_by(df, id) |>
      summarize(admissible = all(admissible))
    inadmissible.ids <- ids.is.admissible[
      !ids.is.admissible$admissible, "id", drop = TRUE
    ]

    df$inadmissible.id <- df$id %in% inadmissible.ids
    E <- mean

  } else {
    df$inadmissible.id <- df$id %in% FALSE
    E <- median
  }

  # ----------------------------------------------------------------------------
  # Bias Plots
  # ----------------------------------------------------------------------------

  dodge <- position_dodge(width = 0.9)
  plot_bias <- function(param = "Y~X:Z", ci.width = 1) {

    tbl <- filter(df,
      !inadmissible.id &
      par %in% param & n == n.i & model.id == model.i & type == type.i
    ) |>
    group_by(
      method, ncat, skew, dist.exo, dist.zeta, par
    ) |>
    summarize(
        bias       = E(bias, na.rm = TRUE),
        se         = sd(est, na.rm = TRUE),
        bias.lower = bias - ci.width * se,
        bias.upper = bias + ci.width * se
    ) |>
    mutate(
      xvar = as.factor(.data[[xv]]),
      par  = sapply(par, \(p) par2tex[[p]])
    )

    plot <- ggplot(tbl, aes(
      x = xvar,
      y = bias,
      colour = method,
      ymin = bias.lower,
      ymax = bias.upper,
      fill = method
    )) +
    geom_col(alpha = 0.2, position = dodge) +
    geom_errorbar(position = dodge, width = 0.25) +
    facet_grid(
      rows = rf,
      cols = cf,
      scales = "fixed",
      labeller = facet_labeller
    ) +
    # ggtitle(sprintf("n = %i, model = %i, %s", n.i, model.i, type.i)) +
    ylab("Bias") +
    xlab(xl) +
    theme_bw()

    min.y <- -0.5
    max.y <- 0.25
    # `na.rm`: sd(est) is NA for any cell left with a single replicate (e.g.
    # when almost everything was inadmissible), which would otherwise make the
    # condition NA rather than FALSE.
    if (any(tbl$bias.lower < min.y, na.rm = TRUE) ||
        any(tbl$bias.upper > max.y, na.rm = TRUE))
      plot <- plot + coord_cartesian(ylim = c(min.y, max.y))

    plot
  }


  # ----------------------------------------------------------------------------
  # SD/SE ratio plots
  # ----------------------------------------------------------------------------
  plot_se_sd_ratio <- function(param = "Y~X:Z") {

    filter(
      df,
      !inadmissible.id &
      par %in% param & n == n.i & model.id == model.i & type == type.i
    ) |>
      group_by(method, ncat, skew, dist.exo, dist.zeta, par) |>
      summarize(
        se = mean(se[admissible.se], na.rm = TRUE),
        sd = sd(est, na.rm = TRUE),
        ratio = se / sd,
        .groups = "drop"
      ) |>
      mutate(
        xvar = as.factor(.data[[xv]]),
        par = sapply(par, \(p) par2tex[[p]])
      ) |>
      ggplot(aes(
        x = xvar,
        y = ratio,
        colour = method,
        shape = method,
        group = method
      )) +
      geom_line(linewidth = 0.5) +
      geom_point(size = 2) +
      facet_grid(
        rows = rf,
        cols = cf,
        scales = "fixed",
        labeller = facet_labeller
      ) +
      ylim(0.8, 1.6) +
      annotate("rect",
        xmin = -Inf, xmax = Inf, ymin = 0.9, ymax = 1.1,
        fill = "grey", alpha = 0.4
      ) +
      # ggtitle(sprintf("n = %i, model = %i, %s", n.i, model.i, type.i)) +
      ylab("SE/SD") +
      xlab(xl) +
      theme_bw()
  }

  # ----------------------------------------------------------------------------
  # SE + SD plots
  # ----------------------------------------------------------------------------
  plot_se_sd <- function(param = "Y~X:Z") {

    filter(
      df,
      !inadmissible.id &
      par == param[[1]] & n == n.i & model.id == model.i & type == type.i
    ) |>
      group_by(method, ncat, skew, dist.exo, dist.zeta, par) |>
      summarize(
        se = mean(se[admissible.se], na.rm = TRUE),
        sd = sd(est, na.rm = TRUE),
        .groups = "drop"
      ) |>
      pivot_longer(
        cols = c("sd", "se"),
        names_to = "measure",
        values_to = "values"
      ) |>
      mutate(
        xvar = as.factor(.data[[xv]]),
        par = sapply(par, \(p) par2tex[[p]])
      ) |>
      ggplot(aes(
        x = xvar,
        y = values,
        colour = method,
        shape = method,
        group = interaction(method, measure),
        linetype = measure
      )) +
      geom_line(linewidth = 0.5, position = position_dodge(width = 0.15)) +
      geom_point(size = 2) +
      facet_grid(
        rows = rf,
        cols = cf,
        scales = "fixed",
        labeller = facet_labeller
      ) +
      # ggtitle(sprintf("n = %i, model = %i, %s", n.i, model.i, type.i)) +
      ylab("SE/SD") +
      xlab(xl) +
      theme_bw()
  }

  # ----------------------------------------------------------------------------
  # Computation Time
  # ----------------------------------------------------------------------------

  dodge <- position_dodge(width = 0.9)
  timeplot <-
    filter(df,
      !inadmissible.id &
      n == n.i & model.id == model.i & type == type.i &
      run.id == timing.runid
    ) |>
    group_by(method, ncat, skew, dist.exo, dist.zeta) |>
    summarize(mean_time = mean(time, na.rm = TRUE)) |>
    mutate(xvar = as.factor(.data[[xv]])) |>
    ggplot(aes(
      x = xvar,
      y = mean_time,
      colour = method,
      fill = method
    )) +
    geom_col(alpha = 0.2, position = dodge) +
    facet_grid(
      rows = rf,
      cols = pf,
      scales = "fixed"
    ) +
    # ggtitle(sprintf("n = %i, model = %i, %s", n.i, model.i, type.i)) +
    ylab("Mean Computation Time (seconds)") +
    xlab(xl) +
    theme_bw()

  # ----------------------------------------------------------------------------
  # Save
  # ----------------------------------------------------------------------------

  plots_nnorm_time[[i]] <- timeplot
  plots_nnorm_bias_l1_l2[[i]] <- plot_bias(c("Y=~y1", "Y=~y2"))
  plots_nnorm_bias_b1[[i]] <- plot_bias("Y~X")
  plots_nnorm_bias_b2[[i]] <- plot_bias("Y~Z")
  plots_nnorm_bias_b1_b2[[i]] <- plot_bias(c("Y~X", "Y~Z"))
  plots_nnorm_bias_b3[[i]] <- plot_bias("Y~X:Z")
  plots_nnorm_se_sd_ratio_b1_b2[[i]] <- plot_se_sd_ratio(c("Y~X", "Y~Z"))
  plots_nnorm_se_sd_ratio_b2[[i]] <- plot_se_sd_ratio("Y~Z")
  plots_nnorm_se_sd_ratio_b3[[i]] <- plot_se_sd_ratio("Y~X:Z")
  plots_nnorm_se_sd_b1[[i]] <- plot_se_sd("Y~X")
  plots_nnorm_se_sd_b2[[i]] <- plot_se_sd("Y~Z")
  plots_nnorm_se_sd_b3[[i]] <- plot_se_sd("Y~X:Z")
  plots_nnorm_inadmissible[[i]] <- pinadmissible
})

dodge <- 0.25
# A plot with computation time accross all conditions and methods. The ordinal
# arm varies with the number of categories; the continuous arm is a single point
# per method, drawn at ncat = 1 so both fit on one axis.
plot_nnorm_time_all <- print(
  df |> mutate(n = as.factor(n), ncat = ifelse(is.na(ncat), 1L, ncat)) |>
  group_by(method, ncat, n) |>
  summarize(mean = mean(time), sd = sd(time), lower = mean - sd, upper = mean + sd) |>
  mutate(ncat = ncat + dodge * as.integer(n=="1000")) |>
  ggplot(aes(x = ncat, y = mean, colour = method, linetype = n, shape = n)) +
  geom_line() + geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper)) + scale_y_log10() +
  theme_bw() +
  ylab("Mean Computation Time (seconds)") +
  xlab("Categories (1 = continuous)")
)

target.type <- "ordinal"
target.n <- 300
target.id <- 1 # currently we only have 1 model in our simulation anyways
idx <- which(simsplit$n == target.n & simsplit$model.id == target.id &
             simsplit$type == target.type)

if (FALSE) {
  print(plots_nnorm_inadmissible[[idx]])
  print(plots_nnorm_time[[idx]])
  print(plots_nnorm_bias_l1_l2[[idx]])
  print(plots_nnorm_bias_b1[[idx]])
  print(plots_nnorm_bias_b2[[idx]])
  print(plots_nnorm_bias_b3[[idx]])
  print(plots_nnorm_se_sd_ratio_b1[[idx]])
  print(plots_nnorm_se_sd_ratio_b2[[idx]])
  print(plots_nnorm_se_sd_ratio_b3[[idx]])
  print(plots_nnorm_se_sd_b1[[idx]])
  print(plots_nnorm_se_sd_b2[[idx]])
  print(plots_nnorm_se_sd_b3[[idx]])
}
