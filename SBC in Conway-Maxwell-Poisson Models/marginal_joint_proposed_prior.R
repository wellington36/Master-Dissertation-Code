library(ggplot2)
library(ggpubr)

# Parameters
N <- 10000 # Number of samples
mu_alpha_prior <- 0.1
mu_beta_prior <- 0.01
nu_alpha_prior <- 1
nu_beta_prior <- 1

# Simulation
set.seed(42) # for reproducibility
mu_samples <- rgamma(N, shape = mu_alpha_prior, rate = mu_beta_prior)
nu_samples <- rgamma(N, shape = nu_alpha_prior, rate = nu_beta_prior)

# Create a data frame
prior_data <- data.frame(mu = mu_samples, nu = nu_samples)

# --- Plotting ---

# 1. Joint Distribution Plot (Scatter Plot)
p_joint <- ggplot(prior_data, aes(x = mu, y = nu)) +
  geom_point(alpha = 0.3, size = 0.5, color = "blue") +
  labs(
    x = expression(paste("Rate Parameter ", mu)),
    y = expression(paste("Dispersion Parameter ", nu)),
    title = "Joint Distribution of Prior Parameters"
  ) +
  theme_minimal() +
  # Limit the x-axis to a reasonable range for visualization
  # The mean is 10, but the variance is 1000, so we'll cap it
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 5)) +
  theme(plot.title = element_text(hjust = 0.5))

# 2. Marginal Distribution for mu (Density Plot)
p_mu_marginal <- ggplot(prior_data, aes(x = mu)) +
  geom_density(fill = "lightblue", alpha = 0.7) +
  labs(x = expression(mu), y = "Density") +
  theme_minimal() +
  coord_cartesian(xlim = c(0, 100)) +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank())

# 3. Marginal Distribution for nu (Density Plot)
p_nu_marginal <- ggplot(prior_data, aes(x = nu)) +
  geom_density(fill = "lightcoral", alpha = 0.7) +
  labs(x = expression(nu), y = "Density") +
  theme_minimal() +
  coord_cartesian(xlim = c(0, 5)) +
  # Flip coordinates for the right margin
  coord_flip() +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())

# Arrange the plots using ggarrange
# The main plot is the joint, with marginals on the top and right
final_plot <- ggarrange(
  p_mu_marginal,
  NULL,
  p_joint,
  p_nu_marginal,
  ncol = 2,
  nrow = 2,
  widths = c(4, 1),
  heights = c(1, 4),
  align = "hv"
)

# Add a common title and save the plot
final_plot_with_title <- annotate_figure(
  final_plot,
  top = text_grob(
    "Marginal and Joint Distribution of the Proposed Prior",
    face = "bold",
    size = 14
  )
)

final_plot_with_title
