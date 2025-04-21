################################################################################
# script to extract raw data and and save forecast data for the case study

################################################################################
## set up

set.seed(91)

library(ncdf4)
library(zoo)
library(crch)
library(ggplot2)
library(evd)
library(pracma)

# source utility functions to load, plot and post-process the raw data
source(here::here("scripts", "case_study_data_utils.R"))


################################################################################
## load data

# load data, remove stations with more than 10% missing values
path <- here::here("data", "tp6_station_")
load_data(path, na_prop = 10)

# plot example of the ensemble forecast trajectories at a random time and station
plot_example(test = F)

# plot map of stations with their corresponding altitudes
plot_map(lonlatalt$lon, lonlatalt$lat, lonlatalt$alt, ymin = -5, ymax = 1600)

# restrict attention to forecasts issued 24h in advance
lead <- 24
lt <- which(lead_times == lead)


################################################################################
## IFS ensemble

# store the IFS ensemble data
fc_ifs <- list(dat = test_fc[, lt, , ])


################################################################################
## smoothed IFS ensemble

# calculate IFS ensemble mean and standard deviation
ens_mn <- apply(fc_ifs$dat, c(1, 2), mean)
ens_sd <- apply(fc_ifs$dat, c(1, 2), sd)

# augment standard deviation when all ensemble members are zero
ens_sd[ens_sd == 0] <- 0.1

# cdf of the logistic distribution censored below at zero
F_x <- function(q, location, scale) pclogis(q, location, scale, left = 0)

# store the forecast cdf with location and scale parameters given by the ensemble mean and standard deviation
fc_smooth <- list(F_x = F_x, location = ens_mn, scale = ens_sd)


################################################################################
## post-processing


### censored logistic

# estimate censored logistic emos parameters at each station (~2 mins)
emos_preds_cl <- lapply(seq_along(stat_ids), function(j) {
  # estimate parameters separately at each station
  st <- stat_ids[j]
  print(paste0('Forecast at Station: ', st, ' (', j, ' from ', length(stat_ids), ')'))

  # extract training data at the station and lead time of interest
  trai <- get_train_data(train_obs, train_fc, train_times, train_years, stat_ids, lead_times[lt], st)

  # extract test data at the station and lead time of interest
  test <- get_test_data(test_obs, test_fc, test_times, stat_ids, lead_times[lt], st)

  # return estimate emos parameters
  fit_emos_crch(trai, test)
})
emos_preds_cl <- simplify2array(emos_preds_cl)

# store the logistic forecast cdf with location and scale parameters given by the estimated emos parameters
fc_emos_cl <- list(F_x = F_x,
                   location = t(emos_preds_cl[, 1, ]),
                   scale = t(emos_preds_cl[, 2, ]))


### censored GEV

# estimate censored GEV emos parameters at each station (~2 hours)
emos_preds_cgev <- lapply(seq_along(stat_ids), function(j) {
  st <- stat_ids[j]
  print(paste0('Forecast at Station: ', st, ' (', j, ' from ', length(stat_ids), ')'))
  trai <- get_train_data(train_obs, train_fc, train_times, train_years, stat_ids, lead_times[lt], st)
  test <- get_test_data(test_obs, test_fc, test_times, stat_ids, lead_times[lt], st)
  fit_emos_cgev(trai, test)
})
emos_preds_cgev <- simplify2array(emos_preds_cgev)

# cdf of the GEV distribution censored below at zero
F_x <- function(q, location, scale, shape) {
  if (length(q) == 1) q <- rep(q, length(location))
  out <- sapply(seq_along(location), function(i) pgev(q[i], location[i], scale[i], shape[i]))
  out[q < 0] <- 0
  return(out)
}

# store the GEV forecast cdf with location, scale and shape parameters given by the estimated emos parameters
fc_emos_cgev <- list(F_x = F_x,
                     location = t(emos_preds_cgev[, 1, ]),
                     scale = t(emos_preds_cgev[, 2, ]),
                     shape = t(emos_preds_cgev[, 3, ]))


################################################################################
## save data

# store auxiliary data (observations, station data, forecast dates, lead times)
aux_data <- list(obs = test_obs[, lt, ],
                 stat_id = stat_ids,
                 stat_info = lonlatalt,
                 time = test_times,
                 lead = lead_times[lt])

# store all forecast (and auxiliary) data in a list
fc_dat <- list(ifs = fc_ifs, smooth = fc_smooth, emos_cl = fc_emos_cl, emos_cgev = fc_emos_cgev, aux_data = aux_data)

# save in an RDS file
saveRDS(fc_dat, file = "scripts/cs_fc_data.RDS")

