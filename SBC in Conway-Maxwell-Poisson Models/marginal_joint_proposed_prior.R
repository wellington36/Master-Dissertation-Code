library(ggplot2)
library(ggpubr)
library(dplyr)

# Parameters
N <- 10000
mu_alpha_prior <- 0.1
mu_beta_prior  <- 0.01
nu_alpha_prior <- 1
nu_beta_prior  <- 1

set.seed(42)

mu_samples <- rgamma(N, shape = mu_alpha_prior, rate = mu_beta_prior)
nu_samples <- rgamma(N, shape = nu_alpha_prior, rate = nu_beta_prior)

prior_data <- data.frame(mu = mu_samples, nu = nu_samples)

# Color palette (match SBC figures)
col_light <- "#CFE8F3"
col_dark  <- "#08519C"

# Joint distribution
p_joint <- ggplot(prior_data, aes(x = mu, y = nu)) +
  geom_point(alpha = 0.15, size = 0.4, color = col_dark) +
  geom_density_2d(color = col_dark, linewidth = 0.6) +
  labs(
    x = expression(paste("Rate Parameter ", mu)),
    y = expression(paste("Dispersion Parameter ", nu))
  ) +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 5)) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

# Marginal for mu
p_mu_marginal <- ggplot(prior_data, aes(x = mu)) +
  geom_density(fill = col_light, color = col_dark, linewidth = 1) +
  labs(x = NULL, y = NULL) +
  coord_cartesian(xlim = c(0, 100)) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.y = element_blank(),
    panel.grid.minor = element_blank()
  )

# Marginal for nu
p_nu_marginal <- ggplot(prior_data, aes(x = nu)) +
  geom_density(fill = col_light, color = col_dark, linewidth = 1) +
  coord_flip(xlim = c(0, 5)) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_blank(),
    panel.grid.minor = element_blank()
  )

# Arrange Layout
final_plot <- ggarrange(
  p_mu_marginal,
  NULL,
  p_joint,
  p_nu_marginal,
  ncol = 2,
  nrow = 2,
  widths = c(4, 1.2),
  heights = c(1.2, 4),
  align = "hv"
)

final_plot_with_title <- annotate_figure(
  final_plot,
  top = text_grob(
    "Marginal and Joint Distribution of the Proposed Prior",
    face = "bold",
    size = 14
  )
)

final_plot_with_title
