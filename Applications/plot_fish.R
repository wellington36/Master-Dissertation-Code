library(brms)

source("rcomp.R")

zinb <- read.csv("https://paul-buerkner.github.io/data/fish.csv")

tab <- table(zinb$count)

x_obs <- as.numeric(names(tab))
y_obs <- as.numeric(tab)
N     <- sum(y_obs)

# ---- Estimated Parameters (Using brms 30k where 20k was warmup) ----
mu_hat     <- 0.001483061
nu_hat     <- 0.08170763
lambda_hat <- mu_hat^nu_hat

lambda_poi_hat <- 1.190736

# COM-Poisson
sim_cmp <- rcomp(
  n      = N,
  lambda = lambda_hat,
  nu     = nu_hat
)

tab_cmp <- table(factor(sim_cmp, levels = x_obs))
y_cmp   <- as.numeric(tab_cmp)

# Poisson
sim_poi <- rpois(
  n      = N,
  lambda = lambda_poi_hat
)

tab_poi <- table(factor(sim_poi, levels = x_obs))
y_poi   <- as.numeric(tab_poi)


# Errors (Observed - Fitted)
err_cmp <- y_obs - y_cmp
err_poi <- y_obs - y_poi

par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))

# ---- Panel 1: Frequencies ----
plot(
  x_obs,
  y_obs,
  type = "b",
  lwd  = 2,
  col  = "steelblue",
  xlab = "Unwanted Pursuit Behaviour (UPB)",
  ylab = "Frequency",
  main = "Observed vs Fitted"
)

lines(x_obs, y_cmp, type = "b", lwd = 2, col = "firebrick")
lines(x_obs, y_poi, type = "b", lwd = 2, col = "darkgreen")

legend(
  "topright",
  legend = c("Observed", "COM-Poisson", "Poisson"),
  col    = c("steelblue", "firebrick", "darkgreen"),
  lwd    = 2,
  bty    = "n"
)

# ---- Panel 2: Errors ----
plot(
  x_obs,
  err_cmp,
  type = "b",
  lwd  = 2,
  col  = "firebrick",
  xlab = "UPB",
  ylab = "Error (Observed - Fitted)",
  main = "Fitting Error",
  ylim = range(c(err_cmp, err_poi))
)

abline(h = 0, lty = 2, col = "gray50")

lines(x_obs, err_poi, type = "b", lwd = 2, col = "darkgreen")

legend(
  "topright",
  legend = c("COM-Poisson Error", "Poisson Error"),
  col    = c("firebrick", "darkgreen"),
  lwd    = 2,
  bty    = "n"
)

# Reset layout
par(mfrow = c(1, 1))
