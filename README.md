# TailCalibration_rep

This repository contains R code to reproduce the results in the following pre-print:

> Allen, S., Koh, J., Segers, J. and Ziegel, J. (2024). 
> Tail calibration of probabilistic forecasts.
> ArXiv preprint.
> [arxiv.org/abs/2407.03167](https://arxiv.org/abs/2407.03167)

The methods introduced in this paper can more generally be implemented using the [`TailCalibration`](https://github.com/sallen12/TailCalibration/tree/main) package in `R`. 

The file `session_info.txt` contains the versions of the packages that were used to generate these results. This was produced using `writeLines(capture.output(sessionInfo()), "session_info.txt")` after sourcing all `R` files in the `scripts` folder.

### Examples

Sourcing the `scripts/examples.R` file generates the plots comprising Figures 3, 4 and 7 of the main paper, and saves them in the `data` directory.


### Simulation studies

Sourcing the `scripts/sim_study_exp.R` file generates the plots comprising Figures 5 and 6 in the main paper, as well as Figure 2 of the supplementary material.

Sourcing the `scripts/sim_study_norm.R` file generates the plots comprising Figure 1 of the supplementary material.

The example and simulation study scripts all run within ~2 minutes.

### Case study

The data employed in the case study (Section 5 of the paper) is publicly available. The `scripts/get_precip_data.py` file can be used to extract the precipitation data, and save it as four `netcdf` files in the `data` directory. This is divided into:
- `tp6_station_1718_fc.ncdf4`: Precipitation ensemble forecasts during the test period of the case study (2017-2018).
- `tp6_station_1718_obs.ncdf4`: The corresponding precipitation observations/measurements during the test period.
- `tp6_station_refo_fc.ncdf4`: Precipitation ensemble re-forecasts, used as training data for the statistical post-processing methods.
- `tp6_station_refo_obs.ncdf4`: The corresponding precipitation observations for the training data set.

Further details about the case study data can be found [here](https://github.com/EUPP-benchmark/climetlab-eumetnet-postprocessing-benchmark). Downloading these four files takes ~30 minutes on a standard (Intel Core i7) CPU. 

The file `scripts/case_study_data.R` pre-processes the data, trains the statistical post-processing models that are compared, and saves the forecasts for the different methods in the file `data/cs_fc_data.RDS`. This pre-processing step takes ~5 minutes, except for the GEV post-processing method, which, without any parallelisation, takes ~2 hours. 

For ease of use, the outputted `data/cs_fc_data.RDS` file is provided in this repository. This file contains a list comprised of the forecast data for the four different forecasting methods compared in the case study (`ifs`, `smooth`, `emos_cl`, `emos_cgev`) and auxiliary data (`aux_data`) such as the observations, station information, and forecast times. The forecasts are either in the form of an ensemble forecast (`dat`), or a forecast distribution function (`F_x`) with additional parameters (e.g. `location`, `scale`, `shape`).

Sourcing the `scripts/case_study_eval.R` file then generates the plots comprising Figure 8 of the main paper, as well as Figure 3 of the supplementary material.

