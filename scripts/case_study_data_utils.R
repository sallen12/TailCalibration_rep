################################################################################
# utility functions to extract raw data and and save forecast data for the case study

################################################################################
## load data

# function to load data
load_data <- function(path, na_prop = 0) {

  ### load training data (20 years of twice-weekly reforecasts, filename = "tp6_station_refo_fc.ncdf4")

  ## reforecast data
  fcst_file <- nc_open(paste0(path, "refo_fc.ncdf4"))
  train_stid_fc <- ncvar_get(fcst_file, varid = "station_id")       # station ids
  train_year_fc <- ncvar_get(fcst_file, varid = "year")             # year of reforecasts
  train_time_fc <- ncvar_get(fcst_file, varid = "time")             # day of the year
  train_lt_fc <- ncvar_get(fcst_file, varid = "step")               # forecast lead time
  train_ens_fc <- ncvar_get(fcst_file, varid = "number")            # ensemble member labels
  train_fc <<- ncvar_get(fcst_file, "tp6")*1000                     # precipitation ensembles (converted from m to mm)
  train_fc <<- aperm(train_fc, c(5, 1, 2, 4, 3))                    # reorder array dimensions
  nc_close(fcst_file)

  ## observation data
  obs_file <- nc_open(paste0(path, "refo_obs.ncdf4"))
  train_stid_obs <- ncvar_get(obs_file, varid = "station_id")       # station ids
  train_year_obs <- ncvar_get(obs_file, varid = "year")             # year of reforecasts
  train_time_obs <- ncvar_get(obs_file, varid = "time")             # day of the year
  train_lt_obs <- ncvar_get(obs_file, varid = "step")               # forecast lead time
  train_obs <<- ncvar_get(obs_file, "tp6")*1000                     # precipitation observations (converted from m to mm)
  lonlatalt <<- data.frame(lon = ncvar_get(obs_file, "longitude"),  # station longitudes
                           lat = ncvar_get(obs_file, "latitude"),   # station latitudes
                           alt = ncvar_get(obs_file, "altitude"),   # station altitudes
                           id = ncvar_get(obs_file, "station_id"))  # station ids
  nc_close(obs_file)



  ### load test data (daily forecasts in 2017 and 2018, filename = "tp6_station_1718_fc.ncdf4")

  ## forecast data
  fcst_file <- nc_open(paste0(path, "1718_fc.ncdf4"))
  test_stid_fc <- ncvar_get(fcst_file, varid = "station_id")        # station ids
  test_time_fc <- ncvar_get(fcst_file, varid = "time") |>           # forecast date
    as.POSIXct(origin = '1970-01-01')                               # convert to date from days since 01-01-1970
  test_lt_fc <- ncvar_get(fcst_file, varid = "step")                # forecast lead time
  test_ens_fc <- ncvar_get(fcst_file, varid = "number")             # ensemble member labels
  test_fc <<- ncvar_get(fcst_file, varid = "tp6")*1000              # precipitation ensembles (converted from m to mm)
  test_fc <<- aperm(test_fc, c(4, 1, 2, 3))                         # reorder array dimensions
  nc_close(fcst_file)

  ## observation data
  obs_file <- nc_open(paste0(path, "1718_obs.ncdf4"))
  test_stid_obs <- ncvar_get(obs_file, varid = "station_id")      # station ids
  test_time_obs <- ncvar_get(obs_file, varid = "time") |>         # forecast date
    as.POSIXct(, origin = '1970-01-01')                           # convert to date from days since 01-01-1970
  test_lt_obs <- ncvar_get(obs_file, varid = "step")              # forecast lead time
  test_obs <<- ncvar_get(obs_file, varid = "tp6")*1000            # precipitation observations (converted from m to mm)
  nc_close(obs_file)


  ### checks

  ## check lead times are the same for forecasts and observations in the training and test data
  if (identical(train_lt_fc, train_lt_obs) &
      identical(test_lt_fc, test_lt_obs) &
      identical(train_lt_fc, test_lt_fc)) {
    lead_times <<- test_lt_obs
  } else {
    stop("Lead times in training and test data do not match")
  }
  ## check station ids are the same for forecasts and observations in the training and test data
  if (identical(train_stid_fc, train_stid_obs) &
      identical(test_stid_fc, test_stid_obs) &
      identical(train_stid_fc, test_stid_fc)) {
    stat_ids <<- test_stid_obs
  } else {
    stop("Station IDs in training and test data do not match")
  }
  ## check forecast times are the same for forecasts and observations
  if (identical(train_time_fc, train_time_obs) &
      identical(test_time_fc, test_time_obs)) {
    train_times <<- train_time_obs
    test_times <<- test_time_obs
  } else {
    stop("Forecast reference times in training and test data do not match")
  }
  ## check reforecast years are the same for forecasts and observations in the training data
  if (identical(train_year_fc, train_year_obs)) {
    train_years <<- train_year_obs
  } else {
    stop("Years in training data do not match")
  }


  train_n_ens <<- length(train_ens_fc)      # ensemble size in the training data
  test_n_ens <<- length(test_ens_fc)        # ensemble size in the test data
  train_obs[is.nan(train_obs)] <<- NA       # assign NA to missing observations in training data
  test_obs[is.nan(test_obs)] <<- NA         # assign NA to missing observations in test data
  train_obs[train_obs < 0] <<- 0            # assign zero to negative precipitation observations in training data
  test_obs[test_obs < 0] <<- 0              # assign zero to negative precipitation observations in test data
  train_fc[train_fc < 0] <<- 0              # assign zero to negative precipitation forecasts in training data
  test_fc[test_fc < 0] <<- 0                # assign zero to negative precipitation forecasts in test data


  ### remove stations with a high proportion of missing observations

  i <- 1
  while (i <= length(stat_ids)) { # cycle through stations
    id <- stat_ids[i]
    train_na <- sapply(1:length(lead_times), function(lt) mean(is.na(train_obs[i, lt, , ])))  # calculate proportion of missing observations at station i in training data
    test_na <- sapply(1:length(lead_times), function(lt) mean(is.na(test_obs[i, lt, ])))      # calculate proportion of missing observations at station i in test data
    if (any(100*train_na > na_prop) || any(100*test_na > na_prop)) {                          # remove stations whose proportion of missing observations is > na_prop
      train_obs <<- train_obs[-i, , , ]
      train_fc <<- train_fc[-i, , , , ]
      test_obs <<- test_obs[-i, , ]
      test_fc <<- test_fc[-i, , , ]
      stat_ids <<- stat_ids[-i]
      lonlatalt <<- lonlatalt[-i, ]
      print(paste("Station", id, "has been removed due to a high proportion of missing values"))
    } else {
      i <- i + 1
    }
  }

}

# function to plot an example ensemble forecast trajectory
plot_example <- function(test = sample(c(T, F), 1)) {
  # sample station at random
  s <- sample(seq_along(stat_ids), 1)
  if (test) {
    n_ens <- test_n_ens
    # sample test date at random
    t <- sample(seq_along(test_times), 1)
    # store ensemble forecast and observation at sampled station and time in a data frame
    df <- data.frame(lt = lead_times,
                     y = c(as.vector(test_fc[s, , t, ]), test_obs[s, , t]),
                     m = as.factor(rep(0:n_ens, each = length(lead_times))))
  } else {
    n_ens <- train_n_ens
    # sample training year and day at random
    y <- sample(seq_along(train_years), 1)
    t <- sample(seq_along(train_times), 1)
    # store ensemble forecast and observation at sampled station and time in a data frame
    df <- data.frame(lt = lead_times,
                     y = c(as.vector(train_fc[s, , y, t, ]), train_obs[s, , y, t]),
                     m = as.factor(rep(0:n_ens, each = length(lead_times))))
  }

  # plot ensemble member and observation trajectories against forecast lead time
  ggplot(df) + geom_line(aes(x = lt, y = y, col = m)) +
    scale_x_continuous(name = "Lead time (hours)", expand = c(0, 0)) +
    scale_y_continuous(name = "Precip. (mm)") +
    scale_color_manual(values = c(rep("grey", n_ens), "black")) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.position = "none")
}

# function to plot values as a function of longitude and latitude
plot_map <- function(lons, lats, z, filename = NULL, ymin = 0, ymax = 15, title = NULL){
  if (is.matrix(z)) z <- as.vector(z)

  # load map data within given longitudes (lons) and latitudes (lats)
  world <- map_data("world")
  ind <- (world$long >= min(lons) & world$long < max(lons)) & (world$lat >= min(lats) & world$lat <= max(lats))
  world <- world[ind, ]

  # plot values of z as a function of longitude and latitude
  df <- data.frame(lat = lats, lon = lons, z = z)
  plot_obj <- ggplot() + geom_point(data = df, aes(lon, lat, fill = z), shape = 21, size = 2) +
    borders("world") +
    coord_fixed(ylim = range(lats), xlim = range(lons)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_gradient(low = "blue", high = "red", limits = c(ymin, ymax),
                        guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
    theme_void() + theme(legend.title = element_blank(), legend.position = "bottom",
                         panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
                         legend.key.width = unit(0.3, "in")) +
    ggtitle(title)


  if (!is.null(filename)) {
    ggsave(filename, plot_obj, width = 2, height = 2.5)
  }

  return(plot_obj)

}


################################################################################
## post-processing

# wrapper to extract training data at a particular lead time (lt) and station (st)
get_train_data <- function(train_obs, train_fc, train_times, train_years, stat_ids, lt, st) {

  trai.all <- NULL

  # cycle over reforecast years
  for (year in train_years) {
    # load observations at chosen station, lead time, and reforecast year
    obs <- train_obs[which(stat_ids == st), which(lead_times == lt), which(train_years == year), ]

    # convert the year and day of the year to a date
    year.dummy <- as.POSIXct(paste0(2017 - 20 + year, '-01-02')) # first Monday in 2017 is 2nd Jan
    init_time <- year.dummy + train_times * 3600 * 24
    time <- init_time + lt * 3600

    # store forecast lead time, station, observations, and dates in a zoo time series
    obs <- cbind(lt, st, obs)
    obs <- zoo(obs, time)
    colnames(obs) <- c('lt', 'stat', 'obs')


    # load IFS ensemble forecast at chosen station, lead time, and reforecast year
    fc <- train_fc[which(stat_ids == st), which(lead_times == lt), which(train_years == year), , ]

    # calculate ensemble mean and standard deviation
    ens.mu <- apply(fc, 1, mean)
    ens.sd <- apply(fc, 1, sd)

    # add forecast parameters to the observation time series
    trai <- cbind(ens.mu, ens.sd)
    trai <- zoo(trai, time)
    trai <- merge(obs, trai)


    # add this time series to the training data from other reforecast years
    trai.all <- rbind(trai.all, trai)
  }

  # remove missing data
  trai <- na.omit(trai.all)

  # restrict training data to data available prior to the beginning of the test data set
  trai <- trai[which(index(trai) < as.POSIXct('2017-01-01')),]

  return(trai)
}

# wrapper to extract test data at a particular lead time (lt) and station (st)
get_test_data <- function(test_obs, test_fc, test_times, stat_ids, lt, st) {

  # load observations at chosen station and lead time
  obs <- test_obs[which(stat_ids == st), which(lead_times == lt), ]

  # store forecast lead time, station, observations, and dates in a zoo time series
  obs <- cbind(lt, st, obs)
  obs <- zoo(obs, test_times)
  colnames(obs) <- c('lt', 'stat', 'obs')

  # load IFS ensemble forecast at chosen station and lead time
  fc <- test_fc[which(stat_ids == st), which(lead_times == lt), , ]

  # calculate ensemble mean and standard deviation
  ens.mu <- apply(fc, 1, mean)
  ens.sd <- apply(fc, 1, sd)

  # add forecast parameters to the observation time series
  test <- cbind(ens.mu, ens.sd)
  test <- zoo(test, test_times)
  test <- merge(obs, test)

  return(test)
}

# wrapper to fit emos models using crch
fit_emos_crch <- function(trai, test, dist = "logistic") {
  # fit parametric distribution (default is censored logistic) to the observations in the training data using crch
  # estimate parameters using the CRPS, revert to maximum likelihood if fitting fails
  fit <- tryCatch(
    {
      fit <- crch(obs ~ ens.mu | ens.sd,
                  data = trai,
                  link.scale = "quadratic",
                  dist = dist,
                  left = 0,
                  type = "crps")

    },
    error = function(cond) {
      message("Fit failed, refitting the model using maximum likelihood")
      fit <- crch(obs ~ ens.mu | ens.sd,
                  data = trai,
                  link.scale = "quadratic",
                  dist = dist,
                  left = 0)
    }
  )
  # return predicted distribution parameters on the test data
  mu <- predict(fit, newdata = test)
  sig <- predict(fit, newdata = test, type = 'scale')
  pred <- cbind(as.numeric(mu), as.numeric(sig))
  colnames(pred) <- c("mu", "sig")
  return(pred)
}

# function to calculate the CRPS of the censored GEV distribution
crps_cgev <- function(y, location, scale, shape, eps = 0.01){
  if (shape >= 1) {
    return(1000) # the CRPS is only valid for shape < 1, assign a high penalty if shape >= 1
  } else if (abs(shape) > eps) {
    # shape != 0

    p0 <- py <- NA*numeric(length(y))
    na_ind <- !is.na(y)
    p0[na_ind] <- pgev(0, loc = location[na_ind], scale = scale[na_ind], shape = shape)
    py[na_ind] <- pgev(y[na_ind], loc = location[na_ind], scale = scale[na_ind], shape = shape)
    p0[p0 < 1e-10] <- 1e-10
    py[py < 1e-10] <- 1e-10

    gamma_1 <- gamma_2 <- NA*numeric(length(y))
    gamma_1[na_ind] <- sapply(py[na_ind], function(x){gammainc(-log(x), 1 - shape)[1]})
    gamma_2[na_ind] <- sapply(p0[na_ind], function(x){gammainc(-2*log(x), 1 - shape)[1]})

    crps1 <- (location - y)*(1 - 2*py)
    crps2 <- location*(p0^2)
    crps3 <- 2*(scale/shape)*(1 - py - gamma_1)
    crps4 <- (scale/shape)*(1 - (p0^2) - (2^shape)*gamma_2)

    crps <- crps1 + crps2 - crps3 + crps4
  } else {
    # shape = 0

    shape_0 <- shape

    shape <- -eps
    p0 <- py <- NA*numeric(length(y))
    na_ind <- !is.na(y)
    p0[na_ind] <- pgev(0, loc = location[na_ind], scale = scale[na_ind], shape = shape)
    py[na_ind] <- pgev(y[na_ind], loc = location[na_ind], scale = scale[na_ind], shape = shape)
    p0[p0 < 1e-10] <- 1e-10
    py[py < 1e-10] <- 1e-10

    gamma_1 <- gamma_2 <- NA*numeric(length(y))
    gamma_1[na_ind] <- sapply(py[na_ind], function(x){gammainc(-log(x), 1 - shape)[1]})
    gamma_2[na_ind] <- sapply(p0[na_ind], function(x){gammainc(-2*log(x), 1 - shape)[1]})

    crps1 <- (location - y)*(1 - 2*py)
    crps2 <- location*(p0^2)
    crps3 <- 2*(scale/shape)*(1 - py - gamma_1)
    crps4 <- (scale/shape)*(1 - (p0^2) - (2^shape)*gamma_2)

    crps_neg <- crps1 + crps2 - crps3 + crps4


    shape <- eps
    p0 <- py <- NA*numeric(length(y))
    na_ind <- !is.na(y)
    p0[na_ind] <- pgev(0, loc = location[na_ind], scale = scale[na_ind], shape = shape)
    py[na_ind] <- pgev(y[na_ind], loc = location[na_ind], scale = scale[na_ind], shape = shape)
    p0[p0 < 1e-10] <- 1e-10
    py[py < 1e-10] <- 1e-10

    gamma_1 <- gamma_2 <- NA*numeric(length(y))
    gamma_1[na_ind] <- sapply(py[na_ind], function(x){gammainc(-log(x), 1 - shape)[1]})
    gamma_2[na_ind] <- sapply(p0[na_ind], function(x){gammainc(-2*log(x), 1 - shape)[1]})

    crps1 <- (location - y)*(1 - 2*py)
    crps2 <- location*(p0^2)
    crps3 <- 2*(scale/shape)*(1 - py - gamma_1)
    crps4 <- (scale/shape)*(1 - (p0^2) - (2^shape)*gamma_2)

    crps_pos <- crps1 + crps2 - crps3 + crps4


    crps <- ((eps - shape_0)*crps_neg + (eps + shape_0)*crps_pos)/(2*eps)
  }
  return(crps)
}

# wrapper to calculate the average CRPS of the censored GEV emos model for given parameter values
crps_gev_obj <- function(par, y, ens.mu, ens.sd) {
  # location parameter is a linear function of the IFS ensemble mean (ens.mu)
  mu <- par[1] + par[2]*ens.mu

  # scale parameter is a linear function of the IFS ensemble standard deviation (ens.sd)
  sig <- (par[3]^2) + (par[4]^2)*ens.sd

  # shape and scale parameters are restricted to be positive using a quadratic link function
  shape <- par[5]^2

  # return average CRPS assigned to the censored GEV forecasts with these parameterss
  crps_cgev(y, location = mu, scale = sig, shape = shape) |> mean()
}

# wrapper to fit emos model with censored GEV distribution
fit_emos_cgev <- function(trai, test) {

  # estimate parameters that minimise the CRPS over the training data
  outpar <- optim(c(1, 1, 1, 1, 0.5), crps_gev_obj, y=trai$obs, ens.mu=trai$ens.mu, ens.sd=trai$ens.sd)$par

  # obtain the predicted location, scale and shape parameters on the test data
  mu <- outpar[1] + outpar[2]*test$ens.mu
  sig <- (outpar[3]^2) + (outpar[4]^2)*test$ens.sd
  shape <- outpar[5]^2

  # return the predicted GEV parameters
  pred <- cbind(as.numeric(mu), as.numeric(sig), as.numeric(shape))
  colnames(pred) <- c("location", "scale", "shape")
  return(pred)
}
