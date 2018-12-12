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
