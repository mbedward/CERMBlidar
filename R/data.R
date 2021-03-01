#' Strata definitions for OFH (Overall Fuel Hazard)
#'
#' @format A data frame (tbl_df) with 4 rows and 3 columns:
#' \describe{
#' \item{name}{Stratum name (NearSurface, Elevated, LowCanopy, UpperCanopy)}
#' \item{lower}{Lower height of stratum (metres)}
#' \item{upper}{Upper height of stratum (metres)}
#' }
"StrataOFH"

#' Strata definitions based on Specht
#'
#' @format A data frame (tbl_df) with 5 rows and 3 columns:
#' \describe{
#' \item{name}{Stratum name (LowShrub, TallShrub, LowTree, MediumTree, TallTree)}
#' \item{lower}{Lower height of stratum (metres)}
#' \item{upper}{Upper height of stratum (metres)}
#' }
"StrataSpecht"

#' Strata definitions for CERMB use.
#'
#' The lowest two strata are 0-0.3m and 0.3-0.5m. Subsequent strata are
#' 0.5m intervals up to 10m; then 1m intervals up to 30m; then 5m
#' intervals up to 80m; then a catch-all layer for >80m.
#'
#' @format A data frame (tbl_df) with 61 rows and 3 columns:
#' \describe{
#' \item{name}{Stratum label ('to0.3', 'to0.5', 'to1.0' ... 'over30.0')}
#' \item{lower}{Lower height of stratum (metres)}
#' \item{upper}{Upper height of stratum (metres)}
#' }
"StrataCERMB"
