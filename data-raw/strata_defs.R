# Create strata definitions and save to package data directory

library(dplyr)

StrataOFH <- tribble(
  ~name,        ~lower, ~upper,
  "Ground",       -Inf,    0.3,
  "NearSurface",   0.3,    0.5,
  "Elevated",      0.5,      5,
  "LowCanopy",       5,     15,
  "UpperCanopy",    15,    Inf
)

usethis::use_data(StrataOFH, overwrite = TRUE)


StrataSpecht <- tribble(
  ~name,        ~lower, ~upper,
  "Ground",       -Inf,    0.5,
  "LowShrub",      0.5,      2,
  "TallShrub",       2,      5,
  "LowTree",         5,     10,
  "MediumTree",     10,     30,
  "TallTree",       30,    Inf
)

usethis::use_data(StrataSpecht, overwrite = TRUE)


# Regular 50cm vertical strata up to 30m

StrataCERMB <- tibble(
  lower = c(-Inf, 0.3, seq(0.5, 30, 0.5)) ) %>%

  mutate(upper = lead(lower),
         name = sprintf("to%.1f", upper)) %>%

  select(name, lower, upper)

i <- nrow(StrataCERMB)
StrataCERMB$name[i] <- sprintf("over%.1f", StrataCERMB$lower[i])
StrataCERMB$upper[i] <- Inf

usethis::use_data(StrataCERMB, overwrite = TRUE)
