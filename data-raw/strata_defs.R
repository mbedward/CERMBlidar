# Create strata definitions and save to package data directory

library(here)
library(dplyr, warn.conflicts = FALSE)


# Overall Fuel Hazard Assessment strata
path <- here("data-raw/strata_ofh.csv")
StrataOFH <- read.csv(path)
usethis::use_data(StrataOFH, overwrite = TRUE)


# Strata adapted from Specht 1970, 1974
path <- here("data-raw/strata_specht.csv")
StrataSpecht <- read.csv(path)
usethis::use_data(StrataSpecht, overwrite = TRUE)


# Original 62-level CERMB strata:
# Regular 50cm layers up to 30m, then a catch-all class
path <- here("data-raw/strata_cermb_old.csv")
StrataCERMB_OLD <- read.csv(path)
usethis::use_data(StrataCERMB_OLD, overwrite = TRUE)


# Updated 52-level CERMB strata:
# Regular 50cm layers up to 10m, then 1m layers up to 30m, then 5m layers up to 80m,
# then a catch-all class
path <- here("data-raw/strata_cermb.csv")
StrataCERMB <- read.csv(path)
usethis::use_data(StrataCERMB, overwrite = TRUE)
