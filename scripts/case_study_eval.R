################################################################################
# script to evaluate the forecasts in the case study

################################################################################
## set up

set.seed(69249837)

library(TailCalibration)
library(crch)
library(evd)
library(ggplot2)

# load forecast data
fc_dat <- readRDS(here::here("data", "cs_fc_data.RDS"))


# remove missing data
y <- fc_dat$aux_data$obs
na_ind <- is.na(y)
y <- y[!na_ind]


# thresholds at which to evaluate tail calibration
y_vec <- seq(0, 20, 0.1)      # thresholds for occurrence ratio diagnostic plots
t_vec <- c(-Inf, 5, 10, 15)   # thresholds for combined and severity ratio diagnostic plots


# specify plot dimensions
width <- 3.5
height <- 3.4


################################################################################
## evaluation


##### IFS

# extract IFS ensemble data and remove forecasts made for missing observations
dat <- fc_dat$ifs$dat
dat <- matrix(dat, ncol = 51)[!na_ind, ]


## combined ratio
com <- tc_prob(y, dat, t = t_vec, lower = 0)
plot_ptc(com, names = c("  S", " 5", "10", "15"), ylims = c(0, 1.9), title = "IFS")
ggsave("plots/cs_com_ifs.png", width = width, height = height)


## severity ratio
sev <- tc_prob(y, dat, t = t_vec, ratio = "sev", lower = 0)
plot_ptc(sev, ratio = "sev", names = c("  S", " 5", "10", "15"), ylims = c(0, 1.02), title = "")
ggsave("plots/cs_sev_ifs.png", width = width, height = height)


## occurrence ratio
occ <- tc_prob(y, dat, t = y_vec, ratio = "occ", lower = 0)
plot_ptc(occ, ratio = "occ", ylims = c(-0.2, 6), xlab = "t (mm)", title = "")
ggsave("plots/cs_occ_ifs.png", width = width, height = height)


## testing

# Kolmogorov-Smirnov test for uniformity of (excess) PIT values
sev <- tc_prob(y, dat, t = t_vec, ratio = "sev", lower = 0, test = TRUE)

# Binomial test for equality of forecast and observed exceedance probabilities
occ <- tc_prob(y, dat, t = t_vec, ratio = "occ", lower = 0, test = TRUE)



##### smoothed IFS

# extract the smoothed IFS forecast distribution, and remove forecasts made for missing observations
F_x <- fc_dat$smooth$F_x                                    # cdf
mu <- as.vector(fc_dat$smooth$location)[!na_ind]            # location parameters
sig <- as.vector(fc_dat$smooth$scale)[!na_ind]              # scale parameters


## combined ratio
com <- tc_prob(y, F_x, t = t_vec, lower = 0, location = mu, scale = sig)
plot_ptc(com, names = c("  S", " 5", "10", "15"), ylims = c(0, 1.9), title = "Smoothed IFS")
ggsave("plots/cs_com_smo.png", width = width, height = height)


## severity ratio
sev <- tc_prob(y, F_x, t = t_vec, ratio = "sev", lower = 0, location = mu, scale = sig)
plot_ptc(sev, ratio = "sev", names = c("  S", " 5", "10", "15"), ylims = c(0, 1.02), title = "")
ggsave("plots/cs_sev_smo.png", width = width, height = height)


## occurrence ratio
occ <- tc_prob(y, F_x, t = y_vec, ratio = "occ", lower = 0, location = mu, scale = sig)
plot_ptc(occ, ratio = "occ", names = c("  S", " 5", "10", "15"), ylims = c(-0.2, 6), xlab = "t (mm)", title = "")
ggsave("plots/cs_occ_smo.png", width = width, height = height)


## testing

# Kolmogorov-Smirnov test for uniformity of (excess) PIT values
sev <- tc_prob(y, F_x, t = t_vec, lower = 0, ratio = "sev", test = TRUE, location = mu, scale = sig)

# Binomial test for equality of forecast and observed exceedance probabilities
occ <- tc_prob(y, F_x, t = t_vec, ratio = "occ", lower = 0, test = TRUE, location = mu, scale = sig)



##### post-processed (logistic)

# extract the censored logistic post-processed forecast distribution, and remove forecasts made for missing observations
F_x <- fc_dat$emos_cl$F_x                                   # cdf
mu <- as.vector(fc_dat$emos_cl$location)[!na_ind]           # location parameters
sig <- as.vector(fc_dat$emos_cl$scale)[!na_ind]             # scale parameters


## combined ratio
com <- tc_prob(y, F_x, t = t_vec, lower = 0, location = mu, scale = sig)
plot_ptc(com, names = c("  S", " 5", "10", "15"), ylims = c(0, 1.9), title = "Post-processed (Logistic)")
ggsave("plots/cs_com_pp.png", width = width, height = height)


## severity ratio
sev <- tc_prob(y, F_x, t = t_vec, lower = 0, ratio = "sev", location = mu, scale = sig)
plot_ptc(sev, ratio = "sev", names = c("  S", " 5", "10", "15"), ylims = c(0, 1.02), title = "")
ggsave("plots/cs_sev_pp.png", width = width, height = height)


## occurrence ratio
occ <- tc_prob(y, F_x, t = y_vec, ratio = "occ", lower = 0, location = mu, scale = sig)
plot_ptc(occ, ratio = "occ", names = c("  S", " 5", "10", "15"), ylims = c(-0.2, 6), xlab = "t (mm)", title = "")
ggsave("plots/cs_occ_pp.png", width = width, height = height)


## testing

# Kolmogorov-Smirnov test for uniformity of (excess) PIT values
sev <- tc_prob(y, F_x, t = t_vec, lower = 0, ratio = "sev", test = TRUE, location = mu, scale = sig)

# Binomial test for equality of forecast and observed exceedance probabilities
occ <- tc_prob(y, F_x, t = t_vec, ratio = "occ", lower = 0, test = TRUE, location = mu, scale = sig)


##### post-processed (GEV)

# extract the censored GEV post-processed forecast distribution, and remove forecasts made for missing observations
F_x <- fc_dat$emos_cgev$F_x                                 # cdf
mu <- as.vector(fc_dat$emos_cgev$location)[!na_ind]         # location parameters
sig <- as.vector(fc_dat$emos_cgev$scale)[!na_ind]           # scale parameters
shape <- as.vector(fc_dat$emos_cgev$shape)[!na_ind]         # shape parameters


## combined ratio
com <- tc_prob(y, F_x, t = t_vec, lower = 0, location = mu, scale = sig, shape = shape)
plot_ptc(com, names = c("  S", " 5", "10", "15"), ylims = c(0, 1.9), title = "Post-processed (GEV)")
ggsave("plots/cs_com_ppgev.png", width = width, height = height)


## severity ratio
sev <- tc_prob(y, F_x, t = t_vec, lower = 0, ratio = "sev", location = mu, scale = sig, shape = shape)
plot_ptc(sev, ratio = "sev", names = c("  S", " 5", "10", "15"), ylims = c(0, 1.02), title = "")
ggsave("plots/cs_sev_ppgev.png", width = width, height = height)


## occurrence ratio
occ <- tc_prob(y, F_x, t = y_vec, ratio = "occ", lower = 0, location = mu, scale = sig, shape = shape)
plot_ptc(occ, ratio = "occ", names = c("  S", " 5", "10", "15"), ylims = c(-0.2, 6), xlab = "t (mm)", title = "")
ggsave("plots/cs_occ_ppgev.png", width = width, height = height)


## testing

# Kolmogorov-Smirnov test for uniformity of (excess) PIT values
sev <- tc_prob(y, F_x, t = t_vec, lower = 0, ratio = "sev", test = TRUE, location = mu, scale = sig, shape = shape)

# Binomial test for equality of forecast and observed exceedance probabilities
occ <- tc_prob(y, F_x, t = t_vec, ratio = "occ", lower = 0, test = TRUE, location = mu, scale = sig, shape = shape)




################################################################################
## evaluation using quantile-based thresholds

# quantile thresholds at which to evaluate tail calibration
q_vec <- c(0.97, 0.99, 0.995)             # thresholds for occurrence ratio diagnostic plots
q_vec_oc <- seq(0.7, 0.999, 0.001)        # thresholds for combined and severity ratio diagnostic plots

# calculate the quantiles of the observations at each quantile threshold
y <- fc_dat$aux_data$obs
a_vec <- sapply(q_vec, function(q) apply(y, 1, quantile, q, na.rm=T) |> rep(ncol(y)))
a_vec_oc <- sapply(q_vec_oc, function(q) apply(y, 1, quantile, q, na.rm=T) |> rep(ncol(y)))

# remove missing observations
y <- y[!na_ind]
a_vec <- a_vec[!na_ind, ]
a_vec_oc <- a_vec_oc[!na_ind, ]

# add the threshold -Inf, which corresponds to standard (S) probabilistic calibration
a_vec <- cbind(-Inf, a_vec)
q_vec <- c(0, q_vec)
q_name <- c(" S", as.character(q_vec))


##### IFS

# extract IFS ensemble data and remove forecasts made for missing observations
dat <- fc_dat$ifs$dat
dat <- matrix(dat, ncol = 51)[!na_ind, ]


## combined ratio
com <- lapply(seq_along(q_vec), function(q) tc_prob(y, dat, t = a_vec[, q], lower = 0, var_t = TRUE))
plot_ptc(com, names = q_name, ylims = c(0, 1.9), title = "IFS")
ggsave("plots/cs_com_ifs_qu.png", width = width, height = height)


## severity ratio
sev <- lapply(seq_along(q_vec), function(q) tc_prob(y, dat, t = a_vec[, q], ratio = "sev", lower = 0, var_t = TRUE))
plot_ptc(sev, ratio = "sev", names = q_name, ylims = c(0, 1.1), title = "")
ggsave("plots/cs_sev_ifs_qu.png", width = width, height = height)


## occurrence ratio
occ <- sapply(seq_along(q_vec_oc), function(q) tc_prob(y, dat, t = a_vec_oc[, q], ratio = "occ", lower = 0, var_t = TRUE))
occ <- data.frame(t = q_vec_oc, rat = occ)
plot_ptc(occ, ratio = "occ", ylims = c(-0.2, 6), xlab = expression(alpha), title = "")
ggsave("plots/cs_occ_ifs_qu.png", width = width, height = height)



##### smoothed IFS

# extract the smoothed IFS forecast distribution, and remove forecasts made for missing observations
F_x <- fc_dat$smooth$F_x                                    # cdf
mu <- as.vector(fc_dat$smooth$location)[!na_ind]            # location parameters
sig <- as.vector(fc_dat$smooth$scale)[!na_ind]              # scale parameters


## combined ratio
com <- lapply(seq_along(q_vec), function(q) tc_prob(y, F_x, t = a_vec[, q], lower = 0, location = mu, scale = sig, var_t = TRUE))
plot_ptc(com, names = q_name, ylims = c(0, 1.9), title = "Smoothed IFS")
ggsave("plots/cs_com_smo_qu.png", width = width, height = height)


## severity ratio
sev <- lapply(seq_along(q_vec), function(q) tc_prob(y, F_x, t = a_vec[, q], ratio = "sev", lower = 0, location = mu, scale = sig, var_t = TRUE))
plot_ptc(sev, ratio = "sev", names = q_name, ylims = c(0, 1.1), title = "")
ggsave("plots/cs_sev_smo_qu.png", width = width, height = height)


## occurrence ratio
occ <- sapply(seq_along(q_vec_oc), function(q) tc_prob(y, F_x, t = a_vec_oc[, q], ratio = "occ", lower = 0, location = mu, scale = sig, var_t = TRUE))
occ <- data.frame(t = q_vec_oc, rat = occ)
plot_ptc(occ, ratio = "occ", ylims = c(-0.2, 6), xlab = expression(alpha), title = "")
ggsave("plots/cs_occ_smo_qu.png", width = width, height = height)



##### post-processed (logistic)

# extract the censored logistic post-processed forecast distribution, and remove forecasts made for missing observations
F_x <- fc_dat$emos_cl$F_x                                   # cdf
mu <- as.vector(fc_dat$emos_cl$location)[!na_ind]           # location parameters
sig <- as.vector(fc_dat$emos_cl$scale)[!na_ind]             # scale parameters


## combined ratio
com <- lapply(seq_along(q_vec), function(q) tc_prob(y, F_x, t = a_vec[, q], lower = 0, location = mu, scale = sig, var_t = TRUE))
plot_ptc(com, names = q_name, ylims = c(0, 1.9), title = "Post-processed (Logistic)")
ggsave("plots/cs_com_pp_qu.png", width = width, height = height)


## severity ratio
sev <- lapply(seq_along(q_vec), function(q) tc_prob(y, F_x, t = a_vec[, q], ratio = "sev", lower = 0, location = mu, scale = sig, var_t = TRUE))
plot_ptc(sev, ratio = "sev", names = q_name, ylims = c(0, 1.1), title = "")
ggsave("plots/cs_sev_pp_qu.png", width = width, height = height)


## occurrence ratio
occ <- sapply(seq_along(q_vec_oc), function(q) tc_prob(y, F_x, t = a_vec_oc[, q], ratio = "occ", lower = 0, location = mu, scale = sig, var_t = TRUE))
occ <- data.frame(t = q_vec_oc, rat = occ)
plot_ptc(occ, ratio = "occ", ylims = c(-0.2, 6), xlab = expression(alpha), title = "")
ggsave("plots/cs_occ_pp_qu.png", width = width, height = height)



##### post-processed (GEV)

# extract the censored GEV post-processed forecast distribution, and remove forecasts made for missing observations
F_x <- fc_dat$emos_cgev$F_x                                 # cdf
mu <- as.vector(fc_dat$emos_cgev$location)[!na_ind]         # location parameters
sig <- as.vector(fc_dat$emos_cgev$scale)[!na_ind]           # scale parameters
shape <- as.vector(fc_dat$emos_cgev$shape)[!na_ind]         # shape parameters


## combined ratio
com <- lapply(seq_along(q_vec), function(q) tc_prob(y, F_x, t = a_vec[, q], lower = 0, location = mu, scale = sig, shape = shape, var_t = TRUE))
plot_ptc(com, names = q_name, ylims = c(0, 1.9), title = "Post-processed (GEV)")
ggsave("plots/cs_com_ppgev_qu.png", width = width, height = height)


## severity ratio
sev <- lapply(seq_along(q_vec), function(q) tc_prob(y, F_x, t = a_vec[, q], ratio = "sev", lower = 0, location = mu, scale = sig, shape = shape, var_t = TRUE))
plot_ptc(sev, ratio = "sev", names = q_name, ylims = c(0, 1.1), title = "")
ggsave("plots/cs_sev_ppgev_qu.png", width = width, height = height)


## occurrence ratio
occ <- sapply(seq_along(q_vec_oc), function(q) tc_prob(y, F_x, t = a_vec_oc[, q], ratio = "occ", lower = 0, location = mu, scale = sig, shape = shape, var_t = TRUE))
occ <- data.frame(t = q_vec_oc, rat = occ)
plot_ptc(occ, ratio = "occ", ylims = c(-0.2, 6), xlab = expression(alpha), title = "")
ggsave("plots/cs_occ_ppgev_qu.png", width = width, height = height)




################################################################################
## normal-based confidence intervals

# function to get confidence intervals for the occurrence ratio
get_occrat_ci <- function(t, y, F_x, ...) {
  exc_F <- 1 - F_x(t, ...)
  exc_y <- y > t
  Fn <- exc_F |> mean()
  Fy <- exc_y |> mean()
  v <- matrix(c(1/Fn, -Fy/(Fn^2)), nrow = 2, ncol = 1)
  Sig <- cov(cbind(exc_y, exc_F))
  sig2 <- (t(v) %*% Sig %*% v) / length(y) |> as.vector()
  return(c(t = t, rat = Fy/Fn, var = sig2))
}

# function to get confidence intervals for the severity ratio
get_sevrat_ci <- function(t, y, F_x, u = seq(0.01, 0.99, 0.01), ...) {
  F_t <- F_x(t, ...)
  Z_F <- (F_x(y, ...) - F_t)/(1 - F_t)
  Z_F[F_t == 1] <- 1
  exc_y <- y > t
  mn_y <- exc_y |> mean()
  n <- length(y)
  if (t < 0) {
    ind <- y == 0
    Z_F[ind] <- runif(sum(ind), 0, Z_F[ind])
  }
  out <- sapply(u, function(u) {
    exc_u <- (Z_F <= u) & exc_y
    mn_u <- exc_u |> mean()
    v <- matrix(c(1/mn_y, -mn_u/(mn_y^2)), nrow = 2, ncol = 1)
    Sig <- cov(cbind(exc_u, exc_y))
    sig2 <- (t(v) %*% Sig %*% v) / n |> as.vector()
    return(c(mean = mn_u/mn_y, var = sig2))
  })
  df <- data.frame(u = u, rat = out[1, ], var = out[2, ])
  return(df)
}

# function to get confidence intervals for the combined ratio
get_comrat_ci <- function(t, y, F_x, u = seq(0.01, 0.99, 0.01), ...) {
  F_t <- F_x(t, ...)
  Z_F <- (F_x(y, ...) - F_t)/(1 - F_t)
  Z_F[F_t == 1] <- 1
  exc_F <- 1 - F_x(t, ...)
  Fn <- exc_F |> mean()
  exc_y <- y > t
  n <- length(y)
  if (t < 0) {
    ind <- y == 0
    Z_F[ind] <- runif(sum(ind), 0, Z_F[ind])
  }
  out <- sapply(u, function(u) {
    exc_u <- (Z_F <= u) & exc_y
    mn_u <- exc_u |> mean()
    v <- matrix(c(1/Fn, -mn_u/(Fn^2)), nrow = 2, ncol = 1)
    Sig <- cov(cbind(exc_u, exc_F))
    sig2 <- (t(v) %*% Sig %*% v) / n |> as.vector()
    return(c(mean = mn_u/Fn, var = sig2))
  })
  df <- data.frame(u = u, rat = out[1, ], var = out[2, ])
  return(df)
}

# function to plot the occurrence ratio with confidence intervals
plot_tc_occ_ci <- function(cal, names = NULL, xlab = NULL, ylab = NULL,
                           xlims = NULL, ylims = NULL, title = NULL, alpha = NULL) {

  if (is.null(xlab)) xlab <- "t"
  if (is.null(ylab)) ylab <- "Occurrence ratio"

  if (is.data.frame(cal)) {
    if (!is.null(alpha)) {
      s <- qnorm((1 + alpha)/2)
      tc <- ggplot(cal) +
        geom_ribbon(aes(x = t, ymin = rat - s*sqrt(var), ymax = rat + s*sqrt(var)), fill = "lightgrey", alpha = 0.5) +
        geom_line(aes(x = t, y = rat))
    } else {
      tc <- ggplot(cal) + geom_line(aes(x = t, y = rat))
    }
  } else {
    df <- do.call(rbind, cal)
    if (!is.null(names)) {
      df$mth <- rep(names, sapply(cal, nrow))
    } else {
      df$mth <- rep(names(cal), sapply(cal, nrow))
    }
    if (!is.null(alpha)) {
      s <- qnorm((1 + alpha)/2)
      tc <- ggplot(df) +
        geom_ribbon(aes(x = t, ymin = rat - s*sqrt(var), ymax = rat + s*sqrt(var), col = mth), alpha = 0.3) +
        geom_line(aes(x = t, y = rat))
    } else {
      tc <- ggplot(df) + geom_line(aes(x = t, y = rat, col = mth))
    }
  }

  tc <- tc + geom_hline(aes(yintercept = 1), linetype = "dotted") +
    scale_x_continuous(name = xlab, limits = xlims, expand = c(0, 0)) +
    scale_y_continuous(name = ylab, limits = ylims, expand = c(0, 0)) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          plot.margin = margin(c(5.5, 10.5, 5.5, 5.5))) +
    ggtitle(title)

  return(tc)
}

# function to plot the combined and severity ratios with confidence intervals
plot_tc_comsev_ci <- function(cal, names = NULL, xlab = NULL, ylab = NULL,
                              xlims = NULL, ylims = NULL, title = NULL, alpha = NULL, com = TRUE) {

  if (is.null(xlab)) xlab <- "u"
  if (is.null(ylab)) if (com) {ylab <- "Combined ratio"} else {ylab = "Severity ratio"}
  if (is.null(xlims)) xlims <- c(0, 1)

  if (is.data.frame(cal)) {
    if (!is.null(alpha)) {
      s <- qnorm((1 + alpha)/2)
      tc <- ggplot(cal) +
        geom_ribbon(aes(x = u, ymin = rat - s*sqrt(var), ymax = rat + s*sqrt(var)), alpha = 0.1) +
        geom_line(aes(x = u, y = rat))
    } else {
      tc <- ggplot(cal) + geom_line(aes(x = u, y = rat))
    }
  } else {
    df <- do.call(rbind, cal)
    if (!is.null(names)) {
      df$mth <- rep(as.factor(names), sapply(cal, nrow))
    } else {
      df$mth <- rep(names(cal), sapply(cal, nrow))
    }
    if (!is.null(alpha)) {
      s <- qnorm((1 + alpha)/2)
      tc <- ggplot(df) +
        geom_ribbon(aes(x = u, ymin = rat - s*sqrt(var), ymax = rat + s*sqrt(var), fill = mth), alpha = 0.1) +
        geom_line(aes(x = u, y = rat, col = mth))
    } else {
      tc <- ggplot(df) + geom_line(aes(x = u, y = rat, col = mth))
    }
  }

  tc <- tc + geom_abline(aes(intercept = 0, slope = 1), linetype = "dotted") +
    scale_x_continuous(name = xlab, limits = xlims, expand = c(0, 0)) +
    scale_y_continuous(name = ylab, limits = ylims, expand = c(0, 0)) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.title = element_blank(),
          legend.justification = c(0, 1),
          legend.position = c(0.01, 0.99),
          plot.margin = margin(c(5.5, 10.5, 5.5, 5.5))) +
    ggtitle(title)

  return(tc)
}

# confidence interval level
alpha <- 0.95


##### IFS

# extract IFS ensemble data and remove forecasts made for missing observations
dat <- fc_dat$ifs$dat
dat <- matrix(dat, ncol = 51)[!na_ind, ]
F_x <- function(x, dat) rowMeans(dat <= x)


## combined ratio
cal <- lapply(t_vec, get_comrat_ci, y = y, F_x = F_x, dat = dat)
names(cal) <- c("  S", " 5", "10", "15")
plot_tc_comsev_ci(cal, alpha = alpha, ylims = c(-0.01, 2.1), title = "IFS")
ggsave("plots/cs_com_ifs_ci.png", width = width, height = height)


## severity ratio
cal <- lapply(t_vec, get_sevrat_ci, y = y, F_x = F_x, dat = dat)
names(cal) <- c("  S", " 5", "10", "15")
plot_tc_comsev_ci(cal, alpha = alpha, ylims = c(-0.01, 1), title = "", com = F)
ggsave("plots/cs_sev_ifs_ci.png", width = width, height = height)


## occurrence ratio
cal <- sapply(y_vec, get_occrat_ci, y = y, F_x = F_x, dat = dat) |> t() |> as.data.frame()
plot_tc_occ_ci(cal, alpha = alpha, ylims = c(-0.2, 7.2), title = "")
ggsave("plots/cs_occ_ifs_ci.png", width = width, height = height)



##### smoothed IFS

# extract the smoothed IFS forecast distribution, and remove forecasts made for missing observations
F_x <- fc_dat$smooth$F_x                                    # cdf
mu <- as.vector(fc_dat$smooth$location)[!na_ind]            # location parameters
sig <- as.vector(fc_dat$smooth$scale)[!na_ind]              # scale parameters


## combined ratio
cal <- lapply(t_vec, get_comrat_ci, y = y, F_x = F_x, location = mu, scale = sig)
names(cal) <- c("  S", " 5", "10", "15")
plot_tc_comsev_ci(cal, alpha = alpha, ylims = c(-0.01, 2.1), title = "Smoothed IFS")
ggsave("plots/cs_com_smo_ci.png", width = width, height = height)


## severity ratio
cal <- lapply(t_vec, get_sevrat_ci, y = y, F_x = F_x, location = mu, scale = sig)
names(cal) <- c("  S", " 5", "10", "15")
plot_tc_sev_ci(cal, alpha = alpha, ylims = c(-0.01, 1), title = "")
ggsave("plots/cs_sev_smo_ci.png", width = width, height = height)


## occurrence ratio
cal <- sapply(y_vec, get_occrat_ci, y = y, F_x = F_x, location = mu, scale = sig) |> t() |> as.data.frame()
plot_tc_occ_ci(cal, alpha = alpha, ylims = c(-0.2, 7.2), title = "")
ggsave("plots/cs_occ_smo_ci.png", width = width, height = height)



##### post-processed (logistic)

# extract the censored logistic post-processed forecast distribution, and remove forecasts made for missing observations
F_x <- fc_dat$emos_cl$F_x                                   # cdf
mu <- as.vector(fc_dat$emos_cl$location)[!na_ind]           # location parameters
sig <- as.vector(fc_dat$emos_cl$scale)[!na_ind]             # scale parameters


## combined ratio
cal <- lapply(t_vec, get_comrat_ci, y = y, F_x = F_x, location = mu, scale = sig)
names(cal) <- c("  S", " 5", "10", "15")
plot_tc_comsev_ci(cal, alpha = alpha, ylims = c(-0.01, 2.1), title = "Post-processed (Logistic)")
ggsave("plots/cs_com_pp_ci.png", width = width, height = height)


## severity ratio
cal <- lapply(t_vec, get_sevrat_ci, y = y, F_x = F_x, location = mu, scale = sig)
names(cal) <- c("  S", " 5", "10", "15")
plot_tc_sev_ci(cal, alpha = alpha, ylims = c(-0.01, 1), title = "")
ggsave("plots/cs_sev_pp_ci.png", width = width, height = height)


## occurrence ratio
cal <- sapply(y_vec, get_occrat_ci, y = y, F_x = F_x, location = mu, scale = sig) |> t() |> as.data.frame()
plot_tc_occ_ci(cal, alpha = alpha, ylims = c(-0.2, 7.2), title = "")
ggsave("plots/cs_occ_pp_ci.png", width = width, height = height)



##### post-processed (GEV)

# extract the censored GEV post-processed forecast distribution, and remove forecasts made for missing observations
F_x <- fc_dat$emos_cgev$F_x                                 # cdf
mu <- as.vector(fc_dat$emos_cgev$location)[!na_ind]         # location parameters
sig <- as.vector(fc_dat$emos_cgev$scale)[!na_ind]           # scale parameters
shape <- as.vector(fc_dat$emos_cgev$shape)[!na_ind]         # shape parameters


## combined ratio
cal <- lapply(t_vec, get_comrat_ci, y = y, F_x = F_x, location = mu, scale = sig, shape = shape)
names(cal) <- c("  S", " 5", "10", "15")
plot_tc_comsev_ci(cal, alpha = alpha, ylims = c(-0.01, 2.1), title = "Post-processed (GEV)")
ggsave("plots/cs_com_ppgev_ci.png", width = width, height = height)


## severity ratio
cal <- lapply(t_vec, get_sevrat_ci, y = y, F_x = F_x, location = mu, scale = sig, shape = shape)
names(cal) <- c("  S", " 5", "10", "15")
plot_tc_sev_ci(cal, alpha = alpha, ylims = c(-0.01, 1), title = "")
ggsave("plots/cs_sev_ppgev_ci.png", width = width, height = height)


## occurrence ratio
cal <- sapply(y_vec, get_occrat_ci, y = y, F_x = F_x, location = mu, scale = sig, shape = shape) |> t() |> as.data.frame()
plot_tc_occ_ci(cal, alpha = alpha, ylims = c(-0.2, 7.2), title = "")
ggsave("plots/cs_occ_ppgev_ci.png", width = width, height = height)


