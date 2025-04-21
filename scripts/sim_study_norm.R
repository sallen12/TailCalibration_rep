################################################################################
# normal simulation study

################################################################################
### set up

library(TailCalibration)
library(ggplot2)


set.seed(631)

# wrapper to save plots
save_plots <- function(prefix, mth = c("cl", "id", "uf", "sr"), width = 3.5, height = 3.4) {
  for (m in mth) {
    # cycle over the three forecasters (climatological, ideal and extremist)
    plotname <- paste0(prefix, "_", m)
    filename <- paste0("plots/sim_norm_", prefix, "_1e6_", m, ".png")
    ggsave(filename, plot = get(plotname), width = width, height = height)
  }
}


################################################################################
### simulate data

N <- 1e6  # sample size

# draw mu values at random from a Norm(0, 1) distribution
mu <- rnorm(N)

# draw observations at random from a Norm(mu, 1) distribution
y <- rnorm(N, mu)

# draw tau values at random from {-1, 1} with probability 0.5
tau <- sample(c(-1, 1), length(mu), replace = T)

# define the cdf of the unfocused forecaster
F_uf <- function(x, m, ta) 0.5*pnorm(x, m) + 0.5*pnorm(x, m + ta)

# quantiles at which tail calibration is to be assessed
rd_q <- c(0.9, 0.95, 0.99, 0.999)
rd_vec <- c(-Inf, quantile(y, rd_q))
names <- c(0, rd_q)

# thresholds at which the occurrence ratio is to be assessed
t_vec <- quantile(y, c(seq(0, 0.99, 0.01), 0.999))


################################################################################
### combined ratio

# get combined ratio values
com_cl <- tc_prob(y, pnorm, t = rd_vec, sd = sqrt(2))                           # climatological forecaster
com_id <- tc_prob(y, pnorm, t = rd_vec, mean = mu)                              # ideal forecaster
com_uf <- tc_prob(y, F_uf, t = rd_vec, m = mu, ta = tau)                        # unfocused forecaster
com_sr <- tc_prob(y, pnorm, t = rd_vec, mean = -mu)                             # sign-reversed forecaster

# plot combined ratio values
com_cl <- com_cl |> plot_ptc(ratio = "com", names = names, ylims = c(0, 1.02), title = "Climatological")
com_id <- com_id |> plot_ptc(ratio = "com", names = names, ylims = c(0, 1.02), title = "Ideal")
com_uf <- com_uf |> plot_ptc(ratio = "com", names = names, ylims = c(0, 1.02), title = "Unfocused")
com_sr <- com_sr |> plot_ptc(ratio = "com", names = names, ylims = c(0, 1.02), title = "Sign-reversed")

# save combined ratio diagnostic plots
save_plots("com")


################################################################################
### severity ratio

# get severity ratio values
sev_cl <- tc_prob(y, pnorm, t = rd_vec, ratio = "sev", sd = sqrt(2))            # climatological forecaster
sev_id <- tc_prob(y, pnorm, t = rd_vec, ratio = "sev", mean = mu)               # ideal forecaster
sev_uf <- tc_prob(y, F_uf, t = rd_vec, ratio = "sev", m = mu, ta = tau)         # unfocused forecaster
sev_sr <- tc_prob(y, pnorm, t = rd_vec, ratio = "sev", mean = -mu)              # sign-reversed forecaster

# plot severity ratio values
sev_cl <- sev_cl |> plot_ptc(ratio = "sev", names = names, ylims = c(0, 1.02), title = "")
sev_id <- sev_id |> plot_ptc(ratio = "sev", names = names, ylims = c(0, 1.02), title = "")
sev_uf <- sev_uf |> plot_ptc(ratio = "sev", names = names, ylims = c(0, 1.02), title = "")
sev_sr <- sev_sr |> plot_ptc(ratio = "sev", names = names, ylims = c(0, 1.02), title = "")

# save severity ratio diagnostic plots
save_plots("sev")


################################################################################
### occurrence ratio

# get occurrence ratio values
occ_cl <- tc_prob(y, pnorm, t = t_vec, ratio = "occ", qu = T, sd = sqrt(2))     # climatological forecaster
occ_id <- tc_prob(y, pnorm, t = t_vec, ratio = "occ", qu = T, mean = mu)        # ideal forecaster
occ_uf <- tc_prob(y, F_uf, t = t_vec, ratio = "occ", qu = T, m = mu, ta = tau)  # unfocused forecaster
occ_sr <- tc_prob(y, pnorm, t = t_vec, ratio = "occ", qu = T, mean = -mu)       # sign-reversed forecaster

# plot occurrence ratio values
occ_cl <- occ_cl |> plot_ptc(ratio = "occ", ylims = c(0, 2), title = "")
occ_id <- occ_id |> plot_ptc(ratio = "occ", ylims = c(0, 2), title = "")
occ_uf <- occ_uf |> plot_ptc(ratio = "occ", ylims = c(0, 2), title = "")
occ_sr <- occ_sr |> plot_ptc(ratio = "occ", ylims = c(0, 2), title = "")

# save occurrence ratio diagnostic plots
save_plots("occ")

