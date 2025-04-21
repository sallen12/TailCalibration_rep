################################################################################
# exponential simulation study


################################################################################
### set up

library(TailCalibration)
library(evd)
library(ggplot2)


set.seed(89)

# wrapper to save plots
save_plots <- function(prefix, mth = c("cl", "id", "ex"), width = 3.5, height = 3.4) {
  for (m in mth) {
    # cycle over the three forecasters (climatological, ideal and extremist)
    plotname <- paste0(prefix, "_", m)
    filename <- paste0("plots/sim_ex_", prefix, "_1e6_", m, ".png")
    ggsave(filename, plot = get(plotname), width = width, height = height)
  }
}


################################################################################
### simulate data

N <- 1e6        # sample size
gamma <- 1/4    # parameter controlling the distribution of delta
v <- 1.4        # parameter controlling the extremist forecaster

# draw delta values at random from a Gamma(1/gamma, 1/gamma) distribution
delta <- rgamma(N, shape = 1/gamma, rate = 1/gamma)

# draw observations at random from an Exponential(delta) distribution
y <- rexp(N, rate = delta)

# quantiles at which tail calibration is to be assessed
rd_q <- c(0.9, 0.95, 0.99, 0.999)
rd_vec <- c(0, quantile(y, rd_q))
names <- c(0, rd_q)

# thresholds at which the occurrence ratio is to be assessed
t_vec <- quantile(y, c(seq(0, 0.99, 0.01), 0.999))


################################################################################
### combined ratio

# get combined ratio values
com_cl <- tc_prob(y, pgpd, t = rd_vec, shape = gamma)                           # climatological forecaster
com_id <- tc_prob(y, pexp, t = rd_vec, rate = delta)                            # ideal forecaster
com_ex <- tc_prob(y, pexp, t = rd_vec, rate = delta/v)                          # extremist forecaster

# plot combined ratio values
com_cl <- com_cl |> plot_ptc(ratio = "com", names = names, ylims = c(0, 1.02), title = "Climatological")
com_id <- com_id |> plot_ptc(ratio = "com", names = names, ylims = c(0, 1.02), title = "Ideal")
com_ex <- com_ex |> plot_ptc(ratio = "com", names = names, ylims = c(0, 1.02), title = "Extremist")

# save combined ratio diagnostic plots
save_plots("com")


################################################################################
### severity ratio

# get severity ratio values
sev_cl <- tc_prob(y, pgpd, t = rd_vec, ratio = "sev", shape = gamma)            # climatological forecaster
sev_id <- tc_prob(y, pexp, t = rd_vec, ratio = "sev", rate = delta)             # ideal forecaster
sev_ex <- tc_prob(y, pexp, t = rd_vec, ratio = "sev", rate = delta/v)           # extremist forecaster

# plot severity ratio values
sev_cl <- sev_cl |> plot_ptc(ratio = "sev", names = names, ylims = c(0, 1.02), title = "")
sev_id <- sev_id |> plot_ptc(ratio = "sev", names = names, ylims = c(0, 1.02), title = "")
sev_ex <- sev_ex |> plot_ptc(ratio = "sev", names = names, ylims = c(0, 1.02), title = "")

# save severity ratio diagnostic plots
save_plots("sev")


################################################################################
### occurrence ratio

# get occurrence ratio values
occ_cl <- tc_prob(y, pgpd, t = t_vec, ratio = "occ", qu = T, shape = gamma)     # climatological forecaster
occ_id <- tc_prob(y, pexp, t = t_vec, ratio = "occ", qu = T, rate = delta)      # ideal forecaster
occ_ex <- tc_prob(y, pexp, t = t_vec, ratio = "occ", qu = T, rate = delta/v)    # extremist forecaster

# plot occurrence ratio values
occ_cl <- occ_cl |> plot_ptc(ratio = "occ", ylims = c(0, 2), title = "")
occ_id <- occ_id |> plot_ptc(ratio = "occ", ylims = c(0, 2), title = "")
occ_ex <- occ_ex |> plot_ptc(ratio = "occ", ylims = c(0, 2), title = "")

# save occurrence ratio diagnostic plots
save_plots("occ")


################################################################################
### conditional combined ratio

# define 3 groups depending on the value of delta
n_grp <- 3
group <- numeric(length(delta))
for (i in 1:n_grp) group[delta >= quantile(delta, (i - 1)/n_grp)] <- paste0("B", i)

# get combined ratio values for each group
com_cl <- tc_cprob(y, pgpd, t = t_vec, group = group, qu = T, shape = gamma)    # climatological forecaster
com_id <- tc_cprob(y, pexp, t = t_vec, group = group, qu = T, rate = delta)     # ideal forecaster
com_ex <- tc_cprob(y, pexp, t = t_vec, group = group, qu = T, rate = delta/v)   # extremist forecaster

# store values for each group in a list
com_b1 <- list("Clim." = com_cl[["B1"]], "Ideal" = com_id[["B1"]], "Extr." = com_ex[["B1"]])
com_b2 <- list("Clim." = com_cl[["B2"]], "Ideal" = com_id[["B2"]], "Extr." = com_ex[["B2"]])
com_b3 <- list("Clim." = com_cl[["B3"]], "Ideal" = com_id[["B3"]], "Extr." = com_ex[["B3"]])

# plot combined ratio values for each group
com_div_b1 <- com_b1 |> plot_tc_sup(ylab = NULL, ylims = c(-0.1, 2.4), title = expression(B[1]))
com_div_b2 <- com_b2 |> plot_tc_sup(ylab = NULL, ylims = c(-0.1, 2.4), title = expression(B[2]))
com_div_b3 <- com_b3 |> plot_tc_sup(ylab = NULL, ylims = c(-0.1, 2.4), title = expression(B[3]))

# save conditional combined ratio plots
save_plots("com_div", mth = c("b1", "b2", "b3"), width = 3, height = 3)

