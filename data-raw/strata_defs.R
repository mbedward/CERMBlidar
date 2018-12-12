# Create strata definitions and save to package data directory

StrataOFH <- dplyr::tribble(
  ~name,        ~lower, ~upper,
  "Ground",       -Inf,    0.3,
  "NearSurface",   0.3,    0.5,
  "Elevated",      0.5,      5,
  "LowCanopy",       5,     15,
  "UpperCanopy",    15,    Inf
)

StrataSpecht <- dplyr::tribble(
  ~name,        ~lower, ~upper,
  "Ground",       -Inf,    0.5,
  "LowShrub",      0.5,      2,
  "TallShrub",       2,      5,
  "LowTree",         5,     10,
  "MediumTree",     10,     30,
  "TallTree",       30,    Inf
)

devtools::use_data(StrataOFH, StrataSpecht, overwrite = TRUE)
