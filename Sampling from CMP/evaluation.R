#set.seed(121)
source("com_poisson_pmf.R")
source("rcomp_exact.R")       # exact sampler
source("rcomp_rou.R")   # RoU sampler
source("rcomp_rejection_benson.R") # rejection sampler

# Dummy sampler for test
rcomp_dummy <- function(n, lambda, nu) {
  rpois(n, lambda)  # just Poisson as placeholder
}

# --- Choose methods to use ---
methods_to_use <- c("rou", "rej_ben", "exa")  
# options: "rou", "rej_ben", "exa", "dummy"

# --- Parameters ---
#lambda <- runif(1, 0.5, 10)
#nu     <- runif(1, 0.5, 2)
lambda <- 2
nu <- 0.5
log_lambda <- log(lambda)
n_sample <- 10000

# --- Storage ---
samples_list <- list()

# --- Simulate conditionally ---
if ("rou" %in% methods_to_use) {
  samples_list$rou <- rcomp_rou(n_sample, lambda, nu)
}
if ("rej_ben" %in% methods_to_use) {
  samples_list$rej_ben <- rcomp_rejection_benson(n_sample, lambda, nu)
}
if ("exa" %in% methods_to_use) {
  samples_list$exa <- rcomp_exact(n_sample, lambda, nu)
}
if ("dummy" %in% methods_to_use) {
  samples_list$dummy <- rcomp_dummy(n_sample, lambda, nu)
}

# --- True pmf (truncate support at reasonable size) ---
K <- max(unlist(samples_list)) + 10
support <- 0:K
log_p <- sapply(support, function(y)
  com_poisson_log_lpmf(y, log_lambda, nu, 1e-12))
p_true <- exp(log_p)
p_true <- p_true / sum(p_true)

df_true <- data.frame(
  x = support,
  y = p_true,        # use 'y' instead of 'prob'
  Method = "True PMF"
)


# --- Empirical histograms ---
make_freq_df <- function(samples, name) {
  tbl <- as.data.frame(table(x = samples))  # explicitly name the column 'x'
  data.frame(
    x = as.numeric(as.character(tbl$x)),
    y = tbl$Freq / sum(tbl$Freq),
    Method = name
  )
}

df_plot <- list(df_true)

if ("rou" %in% methods_to_use) {
  df_plot <- c(df_plot, list(make_freq_df(samples_list$rou, "Rejection RoU")))
}
if ("rej_ben" %in% methods_to_use) {
  df_plot <- c(df_plot, list(make_freq_df(samples_list$rej_ben, "Rejection Benson")))
}
if ("exa" %in% methods_to_use) {
  df_plot <- c(df_plot, list(make_freq_df(samples_list$exa, "Exact")))
}
if ("dummy" %in% methods_to_use) {
  df_plot <- c(df_plot, list(make_freq_df(samples_list$dummy, "Poisson")))
}

df_plot <- do.call(rbind, df_plot)

# --- Plot ---
library(ggplot2)

ggplot(df_plot, aes(x = x, y = y, color = Method, linetype = Method)) +
  geom_point(size = 2) +
  geom_line(size = 0.5) +
  scale_color_manual(values = c(
    "True PMF" = "black",
    "Rejection RoU" = "blue",
    "Rejection Benson" = "orange",
    "Exact" = "green",
    "Poisson" = "red"
  )) +
  scale_linetype_manual(values = c(
    "True PMF" = "solid",
    "Rejection RoU" = "dashed",
    "Rejection Benson" = "dashed",
    "Exact" = "dashed",
    "Poisson" = "dotted"
  )) +
  labs(title = paste0("COM-Poisson λ=", sprintf("%.2f", lambda), 
                      ", ν=", sprintf("%.2f", nu)),
       x = "y", y = "Probability / Relative Frequency") +
  theme_minimal()

library(dplyr)

# Merge with true pmf
df_diff <- df_plot %>%
  left_join(df_true %>% select(x, y_true = y), by = "x") %>%
  mutate(diff = y - y_true) %>%
  filter(Method != "True PMF")

ggplot(df_diff, aes(x = x, y = diff, color = Method)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(size = 2) +
  geom_line() +
  labs(
    title = paste0("Deviation from True COM-Poisson PMF (λ=", sprintf("%.2f", lambda), 
                   ", ν=", sprintf("%.2f", nu), ")"),
    x = "y", y = "Empirical - True"
  ) +
  theme_minimal()


df_ratio <- df_plot %>%
  left_join(df_true %>% select(x, y_true = y), by = "x") %>%
  mutate(ratio = y / y_true) %>%
  filter(Method != "True PMF")

ggplot(df_ratio, aes(x = x, y = ratio, color = Method)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_point(size = 2) +
  geom_line() +
  labs(
    title = paste0("Relative Frequency / True PMF (λ=", sprintf("%.2f", lambda), 
                   ", ν=", sprintf("%.2f", nu), ")"),
    x = "y", y = "Ratio"
  ) +
  theme_minimal()



kl_div <- df_plot %>%
  # join to get y_true from df_true
  left_join(df_true %>% select(x, y_true = y), by = "x") %>%
  filter(Method != "True PMF") %>%
  group_by(Method) %>%
  summarise(
    KL = sum(ifelse(y > 0 & y_true > 0, y * log(y / y_true), 0)),
    .groups = "drop"
  )

ggplot(kl_div, aes(x = Method, y = KL, fill = Method)) +
  geom_col() +
  labs(title = "KL Divergence from True PMF", y = "KL", x = "") +
  theme_minimal()

