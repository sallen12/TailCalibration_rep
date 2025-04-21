################################################################################
###### examples

library(TailCalibration)
library(evd)
library(ggplot2)

set.seed(298301)


################################################################################
###### Example 15: Non-random forecaster

# thresholds and u values at which to assess tail calibration
t <- seq(0, 10, 0.01)
u <- seq(0.01, 0.99, 0.01)

# define a forecast distribution as a deterministic combination of two GPDs

mu2 <- 1                    # location parameter of the second GPD
sig2 <- 0.8                 # scale parameter of the second GPD
th <- mu2 / (1 - sig2)      # threshold at which the GPD densities are equal

# forecast cdf
F_x <- function(x, lower.tail = TRUE) {
  Fcst <- pgpd(x, 0, 1, 1/4)                          # first GPD is implemented at values >= th
  Fcst[x < th] <- pgpd(x[x < th], mu2, sig2, 1/4)     # second GPD is implemented at values < th
  if (lower.tail) {
    return(Fcst)
  } else {
    return(1 - Fcst)
  }
}

# forecast quantile function
F_inv <- function(u) {
  u_th <- pgpd(th, 0, 1, 1/4)
  q <- qgpd(u, 0, 1, 1/4)
  if (any(u < u_th)) {
    q[u < u_th] <- qgpd(u[u < u_th], mu2, sig2, 1/4)
  }
  return(q)
}

# quantile function of the conditional forecast distribution above the threshold t
F_t_inv <- function(u, t) F_inv(u * (1 - F_x(t)) + F_x(t)) - t

# calculate the forecast survival function at each threshold of interest
denom <- F_x(t, lower.tail = FALSE)

# calculate the probability that F_t(Y) <= u and Y > t, where Y follows a standard GPD with shape 1/4
numer <- sapply(u, function(uu) pgpd(t + F_t_inv(uu, t), 0, 1, 1/4) - pgpd(t, 0, 1, 1/4))

# indices of the threshold vector t at which to show results
ind <- c(1, 201, 401, 601)


##### combined ratio

# get combined ratio values and store in a data frame for each threshold
rat <- numer[ind, ]/denom[ind]
cal <- lapply(1:length(ind), function(i) data.frame(u = u, rat = rat[i, ]))
names(cal) <- t[ind]

# plot combined ratio values
plot_ptc(cal)

# save combined ratio diagnostic plot
ggsave("plots/ex_non_com.png", width = 3.2, height = 3)


##### severity ratio

# get severity ratio values and store in a data frame for each threshold
rat <- numer[ind, ] / pgpd(t[ind], 0, 1, 1/4, lower.tail = FALSE)
cal <- lapply(1:length(ind), function(i) data.frame(u = u, rat = rat[i, ]))
names(cal) <- t[ind]

# plot severity ratio values
plot_ptc(cal, ratio = "sev")

# save severity ratio diagnostic plot
ggsave("plots/ex_non_sev.png", width = 3.2, height = 3)


##### occurrence ratio

# get occurrence ratio values and store in a data frame
rat <- pgpd(t, 0, 1, 1/4, lower.tail = FALSE) / denom
df <- data.frame(t = t, r = rat)

# plot occurrence ratio values
plot_ptc(df, ratio = "occ", ylims = c(0, 2))

# save occurrence ratio diagnostic plot
ggsave("plots/ex_non_occ.png", width = 3.2, height = 3)



################################################################################
###### Example 16: Uniform unfocused forecaster

# thresholds and u values at which to assess tail calibration
t <- seq(-1, 0.99, 0.01)
u <- seq(0, 1, 0.01)

# calculate the forecast survival function at each threshold of interest
E_F <- function(y) {
  F_y <- numeric(length(y))
  F_y[y <= -1] <- 0
  F_y[y >= 2] <- 1
  ind <- y > -1 & y <= 1
  F_y[ind] <- (y[ind] + 1)/4
  ind <- y > 1 & y < 2
  F_y[ind] <- (y[ind] + 2)/4
  return(F_y)
}
denom <- 1 - E_F(t)

# calculate the probability that F_t(Y) <= u and Y > t, where Y follows a standard uniform distribution
Q_t <- function(u, t) {
  if (t <= -1) {
    u
  } else if (t > -1 & t <= 0) {
    (pmin(2*u, 1) + pmax(u*(1 - t) + t, 0))/2
  } else if (t > 0) {
    (pmin(u*(2 - t)/(1 - t), 1) + u)/2
  } else {
    NA
  }
}
numer <- t(sapply(t, Q_t, u = u))

# calculate the probability that Y > t, where Y follows a standard uniform distribution
G_t <- pmax(pmin(1 - t, 1), 0)

# indices of the threshold vector t at which to show results
ind <- c(1, 51, 101, 196)


##### combined ratio

# get combined ratio values and store in a data frame for each threshold
rat <- G_t[ind] * numer[ind, ]/denom[ind]
cal <- lapply(1:length(ind), function(i) data.frame(u = u, rat = rat[i, ]))

# plot combined ratio values
plot_ptc(cal, ylims = c(0, 1.4), names = as.factor(t[ind]))

# save combined ratio diagnostic plot
ggsave("plots/ex_uuf_com.png", width = 3.2, height = 3)


##### severity ratio

# get severity ratio values and store in a data frame for each threshold
rat <- numer[ind, ]
cal <- lapply(1:length(ind), function(i) data.frame(u = u, rat = rat[i, ]))

# plot severity ratio values
plot_ptc(cal, ratio = "sev", names = as.factor(t[ind]))

# save severity ratio diagnostic plot
ggsave("plots/ex_uuf_sev.png", width = 3.2, height = 3)


##### occurrence ratio

# get occurrence ratio values and store in a data frame
rat <- G_t / denom
df <- data.frame(t = t, r = rat)

# plot occurrence ratio values
plot_ptc(df, ratio = "occ", ylims = c(-0.2, 2.2))

# save occurrence ratio diagnostic plot
ggsave("plots/ex_uuf_occ.png", width = 3.2, height = 3)



################################################################################
###### Example 18: Optimistic forecaster

N <- 1e6        # sample size
gamma <- 1/4    # parameter controlling the distribution of delta

# draw delta values at random from a Gamma(1/gamma, 1/gamma) distribution
delta <- rgamma(N, shape = 1/gamma, rate = 1/gamma)

# sample X values at random from an Exponential(delta) distribution
x1 <- rexp(N, delta)
x2 <- 2*x1

# sample L values at random from a GPD(1, gamma/2) distribution (independent of delta and X)
L <- rgpd(N, shape = gamma/2)

# define Y as the second largest of (X, 2X, L)
y <- pmax(x1, pmin(x2, L))


# quantiles at which tail calibration is to be assessed
rd_q <- c(0.9, 0.95, 0.99, 0.999)
rd_vec <- c(0, quantile(y, rd_q))
names <- c(0, rd_q)

# thresholds at which the occurrence ratio is to be assessed
t_vec <- quantile(y, c(seq(0, 0.99, 0.01), 0.999))


##### unconditional combined ratio diagnostic plot

# get combined ratio values for an Exponential(delta) forecast distribution
com <- tc_prob(y, pexp, t = rd_vec, qu = T, rate = delta)

# plot combined ratio values
com <- com |> plot_ptc(names = names, ylims = c(0, 1.3), title = "")

# save combined ratio diagnostic plot
ggsave(plot = com, "plots/ex_opt_com_1e6.png", width = 3, height = 3)


##### conditional combined ratio diagnostic plots

# define 3 groups depending on the value of delta
n_grp <- 3
group <- numeric(length(delta))
for (i in 1:n_grp) group[delta >= quantile(delta, (i - 1)/n_grp)] <- paste0("B", i)

# get combined ratio values for an Exponential(delta) forecast distribution for each group
com_div <- tc_cprob(y, pexp, t = t_vec, group = group, qu = T, rate = delta)

# plot and save combined ratio diagnostic plots for each group
com_div[["B1"]] |> plot_tc_sup(ylab = NULL, ylims = c(-0.1, 1.0), title = expression(B[1]))
ggsave("plots/ex_opt_com_div_1e6_B1.png", width = 3, height = 3)
com_div[["B2"]] |> plot_tc_sup(ylab = NULL, ylims = c(-0.1, 1.0), title = expression(B[2]))
ggsave("plots/ex_opt_com_div_1e6_B2.png", width = 3, height = 3)
com_div[["B3"]] |> plot_tc_sup(ylab = NULL, ylims = c(-0.1, 1.0), title = expression(B[3]))
ggsave("plots/ex_opt_com_div_1e6_B3.png", width = 3, height = 3)

