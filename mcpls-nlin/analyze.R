devtools::load_all()
set_project_root()

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(scales)

rdir <- "mcpls-nlin/results/"
files <- dir(rdir)
files <- files[endsWith(files, ".csv")]
paths <- paste0(rdir, files)
ids   <- abbreviate(files, 6)

read <- function(path) {
  name <- last(stringr::str_split_1(path, "/"))
  df <- read.csv(path)
  id <- ids[[name]]
  df$id <- paste0(id, "-", df$id)
  df
}


getModelLoading <- function(par, true) {
  idx <- which.max(par == "X=~x1")
  true[[idx]]
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
    loadings = getModelLoading(par, true),
    admissible.se = all(admissible) & !any(is.na(se) | se > 1) # check SEs when checking admissiblity
  )

# We try to split the simulations into sample size, and model parameter combos
# Within each we look at the performance
simsplit <- expand.grid(
  n = sort(unique(df$n)),
  loadings = sort(unique(df$loadings))
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
admissible <- group_by(df, id, method, loadings, ncat, skew, n) |> 
  summarize(admissible = unique(admissible)) |>
  group_by(method, loadings, ncat, skew, n) |> 
  summarize(nruns = length(admissible),
            ninadmissible = sum(!admissible),
            pinadmissible = sum(!admissible)/length(admissible))

print(admissible, n = 500)


EMPTY_LIST <- vector("list", NROW(simsplit))

plots_inadmissible   <- EMPTY_LIST
plots_time           <- EMPTY_LIST
plots_bias_l1        <- EMPTY_LIST
plots_bias_b1        <- EMPTY_LIST
plots_bias_b2        <- EMPTY_LIST
plots_bias_b3        <- EMPTY_LIST
plots_se_sd_ratio_b1 <- EMPTY_LIST
plots_se_sd_ratio_b2 <- EMPTY_LIST
plots_se_sd_ratio_b3 <- EMPTY_LIST
plots_se_sd_b1       <- EMPTY_LIST
plots_se_sd_b2       <- EMPTY_LIST
plots_se_sd_b3       <- EMPTY_LIST

for (i in seq_len(NROW(simsplit))) suppressMessages({
  cat(sprintf("%i...\n", i))
  # ----------------------------------------------------------------------------
  # Simulation settings
  # ----------------------------------------------------------------------------

  n.i <- simsplit$n[[i]]
  loadings.i <- simsplit$loadings[[i]]

  # ----------------------------------------------------------------------------
  # Inadmissible Solutions
  # ----------------------------------------------------------------------------

  dodge <- position_dodge(width = 0.9)
  pinadmissible <- admissible |>
    filter(n == n.i, loadings == loadings.i) |>
    ggplot(aes(x = ncat, y = pinadmissible, colour = method, fill = method)) +
    geom_col(alpha = 0.2, position = dodge) +
    facet_grid(rows = vars(skew), scales = "fixed") +
    coord_cartesian(ylim = c(0, 1)) +
    scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
    ggtitle(sprintf("Percentage inadmissible solutions (n=%i) (lambda=%.1f)", n.i, loadings.i)) +
    ylab("Percentage inadmissible solutions") +
    xlab("Categories") +
    theme_bw()

  # ----------------------------------------------------------------------------
  # Bias Plots
  # ----------------------------------------------------------------------------

  dodge <- position_dodge(width = 0.9)
  plot_bias <- function(param = "Y~X:Z", ci.width = 1) {
  
    filter(df,
      admissible & par == param[[1]] & n == n.i & loadings == loadings.i
    ) |>
    group_by(
      method, ncat, skew, par
    ) |>
    summarize(
        bias       = mean(bias, na.rm = TRUE),
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
    ggtitle(sprintf("n = %i, loadings = %.1f", n.i, loadings.i)) +
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
      admissible.se & par == param[[1]] & n == n.i & loadings == loadings.i
    ) |>
      group_by(method, ncat, skew, par) |>
      summarize(
        se = mean(se, na.rm = TRUE),
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
      ggtitle(sprintf("n = %i, loadings = %.1f", n.i, loadings.i)) +
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
      admissible.se & par == param[[1]] & n == n.i & loadings == loadings.i
    ) |>
      group_by(method, ncat, skew, par) |>
      summarize(
        se = mean(se, na.rm = TRUE),
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
      ggtitle(sprintf("n = %i, loadings = %.1f", n.i, loadings.i)) +
      ylab("SE/SD") +
      xlab("Categories") +
      theme_bw()
  }

  # ----------------------------------------------------------------------------
  # Computation Time
  # ----------------------------------------------------------------------------

  dodge <- position_dodge(width = 0.9)
  timeplot <-  
    filter(df, admissible & n == n.i & loadings == loadings.i) |>
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
    ggtitle(sprintf("n = %i, loadings = %.1f", n.i, loadings.i)) +
    ylab("Mean Computation Time (seconds)") +
    xlab("Categories") +
    theme_bw()

  # ----------------------------------------------------------------------------
  # Save
  # ----------------------------------------------------------------------------
 
  plots_time[[i]] <- timeplot
  plots_bias_l1[[i]] <- plot_bias("Y=~y1")
  plots_bias_b1[[i]] <- plot_bias("Y~X")
  plots_bias_b2[[i]] <- plot_bias("Y~Z")
  plots_bias_b3[[i]] <- plot_bias("Y~X:Z")
  plots_se_sd_ratio_b1[[i]] <- plot_se_sd_ratio("Y~X")
  plots_se_sd_ratio_b2[[i]] <- plot_se_sd_ratio("Y~Z")
  plots_se_sd_ratio_b3[[i]] <- plot_se_sd_ratio("Y~X:Z")
  plots_se_sd_b1[[i]] <- plot_se_sd("Y~X")
  plots_se_sd_b2[[i]] <- plot_se_sd("Y~Z")
  plots_se_sd_b3[[i]] <- plot_se_sd("Y~X:Z")
  plots_inadmissible[[i]] <- pinadmissible
})

target.n <- 500
target.l <- 0.8
idx <- which(simsplit$n == target.n & simsplit$loadings == target.l)
print(plots_inadmissible[[idx]])
print(plots_time[[idx]])
print(plots_bias_l1[[idx]])
print(plots_bias_b1[[idx]])
print(plots_bias_b2[[idx]])
print(plots_bias_b3[[idx]])
print(plots_se_sd_ratio_b1[[idx]])
print(plots_se_sd_ratio_b2[[idx]])
print(plots_se_sd_ratio_b3[[idx]])
print(plots_se_sd_b1[[idx]])
print(plots_se_sd_b2[[idx]])
print(plots_se_sd_b3[[idx]])
