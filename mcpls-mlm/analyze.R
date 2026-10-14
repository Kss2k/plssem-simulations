devtools::load_all()
set_project_root()

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(scales)

testfiles <- FALSE
rdir <- "mcpls-mlm/results/"
time.run.id <- "v1-tuf" # computation times are only compared within one machine
files <- dir(rdir)
files <- files[endsWith(files, ".csv")]
if (!testfiles) files <- files[!grepl("test", files) & !grepl("extra", files)]

paths <- paste0(rdir, files)
ids   <- stats::setNames(paste(
  stringr::str_extract(files, "^results-(v1-[A-Za-z]+)-([0-9]+)-", group = 1),
  stringr::str_extract(files, "^results-(v1-[A-Za-z]+)-([0-9]+)-", group = 2),
  sep = "-"
), nm = files)

read <- function(path) {
  name <- last(stringr::str_split_1(path, "/"))
  df <- read.csv(path)
  id <- ids[[name]]
  df$id <- paste0(id, "-", df$id)
  df
}


methods_ordered <- c("MC-PLS", "Mplus")

# Only the standardized estimates have true values. The random slope SDs
# (`rsd.*`) are only estimated by MC-PLS (they are always NA for Mplus)
df <- do.call(rbind, lapply(paths, read)) |>
  filter(type == "standardized", !(method == "Mplus" & startsWith(par, "rsd."))) |>
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
  k = sort(unique(df$k)),
  model.id = sort(unique(df$model.id)),
  drop.inadmissible = c(TRUE, FALSE)
)

par2tex <- list(
  `Yw~Xw`     = "gamma[1]^W",
  `Yw~Zw`     = "gamma[2]^W",
  `Yw~Xw:Zw`  = "gamma[3]^W",
  `Xw~~Zw`    = "phi^W",
  `Yb~Xb.l2`  = "gamma[1]^B",
  `Yb~Zb.l2`  = "gamma[2]^B",
  `Xb~~Zb.l2` = "phi^B",
  `s~Zb.l2`   = "gamma[s]^B",
  `rsd.s`     = "sigma[s]",
  `Yw=~y1`    = "lambda[7]^W",
  `Yw=~y2`    = "lambda[8]^W",
  `Yb=~y1.l2` = "lambda[7]^B",
  `Yb=~y2.l2` = "lambda[8]^B"
)

# Count inadmissibles
admissible <- group_by(df, id, method, model.id, ncat, skew, n, k) |> 
  summarize(admissible = unique(admissible)) |>
  group_by(method, model.id, ncat, skew, n, k) |> 
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
plots_bias_between      <- EMPTY_LIST
plots_bias_rslope       <- EMPTY_LIST
plots_se_sd_ratio_b1    <- EMPTY_LIST
plots_se_sd_ratio_b2    <- EMPTY_LIST
plots_se_sd_ratio_b3    <- EMPTY_LIST
plots_se_sd_ratio_b1_b2 <- EMPTY_LIST
plots_se_sd_ratio_between <- EMPTY_LIST
plots_se_sd_ratio_rslope  <- EMPTY_LIST
plots_se_sd_b1          <- EMPTY_LIST
plots_se_sd_b2          <- EMPTY_LIST
plots_se_sd_b3          <- EMPTY_LIST

for (i in seq_len(NROW(simsplit))) suppressMessages({
  cat(sprintf("%i...\n", i))
  # ----------------------------------------------------------------------------
  # Simulation settings
  # ----------------------------------------------------------------------------

  n.i <- simsplit$n[[i]]
  k.i <- simsplit$k[[i]]
  model.i <- simsplit$model.id[[i]]
  drop.inadmissible <- simsplit$drop.inadmissible[[i]]

  # ----------------------------------------------------------------------------
  # Inadmissible Solutions
  # ----------------------------------------------------------------------------

  dodge <- position_dodge(width = 0.9)
  pinadmissible <- admissible |>
    filter(n == n.i, k == k.i, model.id == model.i) |>
    ggplot(aes(x = ncat, y = pinadmissible, colour = method, fill = method)) +
    geom_col(alpha = 0.2, position = dodge) +
    facet_grid(rows = vars(skew), scales = "fixed") +
    # coord_cartesian(ylim = c(0, 1)) +
    scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
    # ggtitle(sprintf("Percentage inadmissible solutions (n=%i) model %d", n.i, model.i)) +
    ylab("Percentage inadmissible solutions") +
    xlab("Categories") +
    theme_bw()
  
  pinadmissible2 <- 
    admissible |> mutate(
      ncat = as.factor(ncat),
      pinadmissible = 100 * pinadmissible,
      pinadmissible.scaled = (pinadmissible - max(pinadmissible))^(1/8)
    ) |>
    filter(n == n.i, k == k.i, model.id == model.i) |>
    ggplot(aes(x = ncat, y = skew)) +
    geom_tile(aes(fill=pinadmissible)) +
    geom_text(aes(label=paste0(round(pinadmissible,1), "%"))) +
    facet_wrap(~method) +
    ylab("Percentage inadmissible solutions") +
    xlab("Categories") +
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
  plot_bias <- function(param = "Yw~Xw:Zw", ci.width = 1) {
  
    tbl <- filter(df,
      !inadmissible.id &
      par %in% param & n == n.i & k == k.i & model.id == model.i 
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
    )

    plot <- ggplot(tbl, aes(
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

    min.y <- -0.5
    max.y <- 0.25
    if (any(tbl$bias.lower < min.y) || any(tbl$bias.upper > max.y))
      plot <- plot + coord_cartesian(ylim = c(min.y, max.y))

    plot
  }


  # ----------------------------------------------------------------------------
  # SD/SE ratio plots
  # ----------------------------------------------------------------------------
  plot_se_sd_ratio <- function(param = "Yw~Xw:Zw") {

    filter(
      df,
      !inadmissible.id &
      par %in% param & n == n.i & k == k.i & model.id == model.i 
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
      # ggtitle(sprintf("n = %i, model = %i", n.i, model.i)) +
      ylab("SE/SD") +
      xlab("Categories") +
      theme_bw()
  }

  # ----------------------------------------------------------------------------
  # SE + SD plots
  # ----------------------------------------------------------------------------
  plot_se_sd <- function(param = "Yw~Xw:Zw") {

    filter(
      df,
      !inadmissible.id &
      par == param[[1]] & n == n.i & k == k.i & model.id == model.i 
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
      # ggtitle(sprintf("n = %i, model = %i", n.i, model.i)) +
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
      n == n.i & k == k.i & model.id == model.i &
      grepl(time.run.id, id)
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
    # ggtitle(sprintf("n = %i, model = %i", n.i, model.i)) +
    ylab("Mean Computation Time (seconds)") +
    xlab("Categories") +
    theme_bw()

  # ----------------------------------------------------------------------------
  # Save
  # ----------------------------------------------------------------------------
 
  plots_time[[i]] <- timeplot
  # Some parameters only exist in one of the models (e.g., `Yw~Xw:Zw` and `rsd.s`)
  has_par <- \(param) any(df$par %in% param & df$model.id == model.i)
  if_par  <- \(param, plot) if (has_par(param)) plot(param) else NULL

  plots_bias_l1_l2[[i]] <- plot_bias(c("Yw=~y1", "Yw=~y2", "Yb=~y1.l2", "Yb=~y2.l2"))
  plots_bias_b1[[i]] <- plot_bias("Yw~Xw")
  plots_bias_b2[[i]] <- plot_bias("Yw~Zw")
  plots_bias_b1_b2[[i]] <- plot_bias(c("Yw~Xw", "Yw~Zw"))
  plots_bias_b3[[i]] <- if_par("Yw~Xw:Zw", plot_bias)
  plots_bias_between[[i]] <- plot_bias(c("Yb~Xb.l2", "Yb~Zb.l2", "Xb~~Zb.l2"))
  plots_bias_rslope[[i]] <- if_par(c("s~Zb.l2", "rsd.s"), plot_bias)
  plots_se_sd_ratio_b1_b2[[i]] <- plot_se_sd_ratio(c("Yw~Xw", "Yw~Zw"))
  plots_se_sd_ratio_b2[[i]] <- plot_se_sd_ratio("Yw~Zw")
  plots_se_sd_ratio_b3[[i]] <- if_par("Yw~Xw:Zw", plot_se_sd_ratio)
  plots_se_sd_ratio_between[[i]] <- plot_se_sd_ratio(c("Yb~Xb.l2", "Yb~Zb.l2", "Xb~~Zb.l2"))
  plots_se_sd_ratio_rslope[[i]] <- if_par(c("s~Zb.l2", "rsd.s"), plot_se_sd_ratio)
  plots_se_sd_b1[[i]] <- plot_se_sd("Yw~Xw")
  plots_se_sd_b2[[i]] <- plot_se_sd("Yw~Zw")
  plots_se_sd_b3[[i]] <- if_par("Yw~Xw:Zw", plot_se_sd)
  plots_inadmissible[[i]] <- pinadmissible
})

dodge <- 0.25
# A plot with computation time accross all conditions and methods
plot_time_all <- print(
  df |> mutate(n = as.factor(n), ncat = factor(ncat)) |> # ncat = Inf (continuous)
  group_by(method, model.id, ncat, n) |>
  summarize(mean = mean(time), sd = sd(time), lower = mean - sd, upper = mean + sd) |>
  ggplot(aes(x = ncat, y = mean, colour = method, linetype = n, shape = n,
             group = interaction(method, n))) +
  geom_line(position = position_dodge(width = dodge)) +
  geom_point(position = position_dodge(width = dodge)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(width = dodge)) +
  scale_y_log10() +
  facet_wrap(~model.id) +
  theme_bw() +
  ylab("Mean Computation Time (seconds)") +
  xlab("Categories")
)

target.n <- 20
target.k <- 30
target.id <- 1 # 1 = within_xz, 2 = random_slope
target.drop <- TRUE # drop inadmissible solutions?
idx <- which(simsplit$n == target.n & simsplit$k == target.k & simsplit$model.id == target.id &
             simsplit$drop.inadmissible == target.drop)

if (FALSE) {
  print(plots_inadmissible[[idx]])
  print(plots_time[[idx]])
  print(plots_bias_l1_l2[[idx]])
  print(plots_bias_b1[[idx]])
  print(plots_bias_b2[[idx]])
  print(plots_bias_b3[[idx]])
  print(plots_bias_between[[idx]])
  print(plots_bias_rslope[[idx]])
  print(plots_se_sd_ratio_b1[[idx]])
  print(plots_se_sd_ratio_b2[[idx]])
  print(plots_se_sd_ratio_b3[[idx]])
  print(plots_se_sd_ratio_between[[idx]])
  print(plots_se_sd_ratio_rslope[[idx]])
  print(plots_se_sd_b1[[idx]])
  print(plots_se_sd_b2[[idx]])
  print(plots_se_sd_b3[[idx]])
}
