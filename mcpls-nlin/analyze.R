devtools::load_all()
set_project_root()

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(scales)

testfiles <- TRUE
rdir <- "mcpls-nlin/results/"
files <- dir(rdir)
files <- files[endsWith(files, ".csv")]
if (!testfiles) files <- files[!grepl("test", files) & !grepl("extra", files)]

paths <- paste0(rdir, files)
ids   <- abbreviate(files, 6)

read <- function(path) {
  name <- last(stringr::str_split_1(path, "/"))
  df <- read.csv(path)
  id <- ids[[name]]
  df$id <- paste0(id, "-", df$id)
  df
}


methods_ordered <- c("PLS", "PLSc", "MC-OrdPLSc", "Mplus")

df <- do.call(rbind, lapply(paths, read)) |>
  mutate(
    bias = est - true,
    method = factor(method,
      levels = methods_ordered,
      labels = methods_ordered
    ),
  ) |>
  group_by(id, method) |>
  mutate(
    parcombo = paste0(paste0(par, "=", true), collapse = ","),
    admissible.se = all(admissible) & !any(is.na(se) | se > 1) # check SEs when checking admissiblity
  )

# We try to split the simulations into sample size, and model parameter combos
# Within each we look at the performance
simsplit <- expand.grid(
  n = sort(unique(df$n)),
  model.id = sort(unique(df$model.id)),
  drop.inadmissible = c(TRUE, FALSE)
)

par2tex <- list(
  `Y~X`   = "gamma[1]",
  `Y~Z`   = "gamma[2]",
  `Y~X:Z` = "gamma[3]",
  `Y=~y1` = "lambda[7]",
  `Y=~y2` = "lambda[8]",
  `Y=~y3` = "lambda[9]"
)

# Count inadmissibles
admissible <- group_by(df, id, method, model.id, ncat, skew, n) |> 
  summarize(admissible = unique(admissible)) |>
  group_by(method, model.id, ncat, skew, n) |> 
  summarize(nruns = length(admissible),
            ninadmissible = sum(!admissible),
            pinadmissible = sum(!admissible)/length(admissible))

print(admissible, n = 500)


EMPTY_LIST <- vector("list", NROW(simsplit))

plots_inadmissible      <- EMPTY_LIST
plots_time              <- EMPTY_LIST
plots_bias_l1_l2        <- EMPTY_LIST
plots_bias_b1           <- EMPTY_LIST
plots_bias_b1_b2        <- EMPTY_LIST
plots_bias_b2           <- EMPTY_LIST
plots_bias_b3           <- EMPTY_LIST
plots_se_sd_ratio_b1    <- EMPTY_LIST
plots_se_sd_ratio_b2    <- EMPTY_LIST
plots_se_sd_ratio_b3    <- EMPTY_LIST
plots_se_sd_ratio_b1_b2 <- EMPTY_LIST
plots_se_sd_b1          <- EMPTY_LIST
plots_se_sd_b2          <- EMPTY_LIST
plots_se_sd_b3          <- EMPTY_LIST

for (i in seq_len(NROW(simsplit))) suppressMessages({
  cat(sprintf("%i...\n", i))
  # ----------------------------------------------------------------------------
  # Simulation settings
  # ----------------------------------------------------------------------------

  n.i <- simsplit$n[[i]]
  model.i <- simsplit$model.id[[i]]
  drop.inadmissible <- simsplit$drop.inadmissible[[i]]

  # ----------------------------------------------------------------------------
  # Inadmissible Solutions
  # ----------------------------------------------------------------------------

  dodge <- position_dodge(width = 0.9)
  pinadmissible <- admissible |>
    filter(n == n.i, model.id == model.i) |>
    ggplot(aes(x = ncat, y = pinadmissible, colour = method, fill = method)) +
    geom_col(alpha = 0.2, position = dodge) +
    facet_grid(rows = vars(skew), scales = "fixed") +
    # coord_cartesian(ylim = c(0, 1)) +
    scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
    ggtitle(sprintf("Percentage inadmissible solutions (n=%i) model %d", n.i, model.i)) +
    ylab("Percentage inadmissible solutions") +
    xlab("Categories") +
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
  
    filter(df,
      !inadmissible.id &
      par %in% param & n == n.i & model.id == model.i 
    ) |>
    group_by(
      method, ncat, skew, par
    ) |>
    summarize(
        bias       = E(bias, na.rm = TRUE),
        se         = sd(est, na.rm = TRUE),
        bias.lower = bias - ci.width * se,
        bias.upper = bias + ci.width * se
    ) |> 
    mutate(
      ncat = as.factor(ncat),
      par  = sapply(par, \(p) par2tex[[p]])
    ) |>
    ggplot(aes(
      x = ncat,
      y = bias,
      colour = method,
      ymin = bias.lower,
      ymax = bias.upper,
      fill = method
    )) +
    geom_col(alpha = 0.2, position = dodge) +
    geom_errorbar(position = dodge, width = 0.25) +
    facet_grid(
      rows = vars(skew),
      cols = vars(par),
      scales = "fixed",
      labeller = label_parsed
    ) +
    # ggtitle(sprintf("n = %i, model = %i", n.i, model.i)) +
    ylab("Bias") +
    xlab("Categories") +
    theme_bw()
  }


  # ----------------------------------------------------------------------------
  # SD/SE ratio plots
  # ----------------------------------------------------------------------------
  plot_se_sd_ratio <- function(param = "Y~X:Z") {

    filter(
      df,
      !inadmissible.id &
      par %in% param & n == n.i & model.id == model.i 
    ) |>
      group_by(method, ncat, skew, par) |>
      summarize(
        se = mean(se[admissible.se], na.rm = TRUE),
        sd = sd(est, na.rm = TRUE),
        ratio = se / sd,
        .groups = "drop"
      ) |>
      mutate(
        ncat = factor(ncat),
        par = sapply(par, \(p) par2tex[[p]])
      ) |>
      ggplot(aes(
        x = ncat,
        y = ratio,
        colour = method,
        shape = method,
        group = method
      )) +
      geom_line(linewidth = 0.5) +
      geom_point(size = 2) +
      facet_grid(
        rows = vars(skew),
        cols = vars(par),
        scales = "fixed",
        labeller = label_parsed
      ) +
      ylim(0.8, 1.6) +
      annotate("rect",
        xmin = -Inf, xmax = Inf, ymin = 0.9, ymax = 1.1, 
        fill = "grey", alpha = 0.4
      ) +
      ggtitle(sprintf("n = %i, model = %i", n.i, model.i)) +
      ylab("SE/SD") +
      xlab("Categories") +
      theme_bw()
  }

  # ----------------------------------------------------------------------------
  # SE + SD plots
  # ----------------------------------------------------------------------------
  plot_se_sd <- function(param = "Y~X:Z") {

    filter(
      df,
      !inadmissible.id &
      par == param[[1]] & n == n.i & model.id == model.i 
    ) |>
      group_by(method, ncat, skew, par) |>
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
        ncat = factor(ncat),
        par = sapply(par, \(p) par2tex[[p]])
      ) |>
      ggplot(aes(
        x = ncat,
        y = values,
        colour = method,
        shape = method,
        group = interaction(method, measure),
        linetype = measure
      )) +
      geom_line(linewidth = 0.5, position = position_dodge(width = 0.15)) +
      geom_point(size = 2) +
      facet_grid(
        rows = vars(skew),
        cols = vars(par),
        scales = "fixed",
        labeller = label_parsed
      ) +
      ggtitle(sprintf("n = %i, model = %i", n.i, model.i)) +
      ylab("SE/SD") +
      xlab("Categories") +
      theme_bw()
  }

  # ----------------------------------------------------------------------------
  # Computation Time
  # ----------------------------------------------------------------------------

  dodge <- position_dodge(width = 0.9)
  timeplot <-  
    filter(df,
      !inadmissible.id &
      n == n.i & model.id == model.i
    ) |>
    group_by(method, ncat, skew) |>
    summarize(mean_time = mean(time, na.rm = TRUE)) |>
    mutate(ncat = as.factor(ncat)) |>
    ggplot(aes(
      x = ncat,
      y = mean_time,
      colour = method,
      fill = method
    )) +
    geom_col(alpha = 0.2, position = dodge) +
    facet_grid(
      rows = vars(skew),
      scales = "fixed"
    ) +
    ggtitle(sprintf("n = %i, model = %i", n.i, model.i)) +
    ylab("Mean Computation Time (seconds)") +
    xlab("Categories") +
    theme_bw()

  # ----------------------------------------------------------------------------
  # Save
  # ----------------------------------------------------------------------------
 
  plots_time[[i]] <- timeplot
  plots_bias_l1_l2[[i]] <- plot_bias(c("Y=~y1", "Y=~y2"))
  plots_bias_b1[[i]] <- plot_bias("Y~X")
  plots_bias_b2[[i]] <- plot_bias("Y~Z")
  plots_bias_b1_b2[[i]] <- plot_bias(c("Y~X", "Y~Z"))
  plots_bias_b3[[i]] <- plot_bias("Y~X:Z")
  plots_se_sd_ratio_b1_b2[[i]] <- plot_se_sd_ratio(c("Y~X", "Y~Z"))
  plots_se_sd_ratio_b2[[i]] <- plot_se_sd_ratio("Y~Z")
  plots_se_sd_ratio_b3[[i]] <- plot_se_sd_ratio("Y~X:Z")
  plots_se_sd_b1[[i]] <- plot_se_sd("Y~X")
  plots_se_sd_b2[[i]] <- plot_se_sd("Y~Z")
  plots_se_sd_b3[[i]] <- plot_se_sd("Y~X:Z")
  plots_inadmissible[[i]] <- pinadmissible
})

target.n <- 300
target.id <- 1 # currently we only have 1 model in our simulation anyways
idx <- which(simsplit$n == target.n & simsplit$model.id == target.id)

if (FALSE) {
  print(plots_inadmissible[[idx]])
  print(plots_time[[idx]])
  print(plots_bias_l1_l2[[idx]])
  print(plots_bias_b1[[idx]])
  print(plots_bias_b2[[idx]])
  print(plots_bias_b3[[idx]])
  print(plots_se_sd_ratio_b1[[idx]])
  print(plots_se_sd_ratio_b2[[idx]])
  print(plots_se_sd_ratio_b3[[idx]])
  print(plots_se_sd_b1[[idx]])
  print(plots_se_sd_b2[[idx]])
  print(plots_se_sd_b3[[idx]])
}
