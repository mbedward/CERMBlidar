#' Strata definitions used for Overall Fuel Hazard Assessment in New South Wales.
#'
#' @format A data frame with 5 rows and 5 columns:
#' \describe{
#' \item{strata_def}{Name for this strata scheme ('ofh')}
#' \item{index}{Integer index, ordered by height level}
#' \item{label}{Stratum name (Ground, NearSurface, Elevated, LowCanopy, UpperCanopy)}
#' \item{lowerht}{Lower height of stratum (metres; -999.0 for ground layer)}
#' \item{upperht}{Upper height of stratum (metres; 999.0 for top layer)}
#' }
"StrataOFH"


#' Strata definitions adapted from Specht 1970, 1974.
#'
#' @format A data frame with 6 rows and 5 columns:
#' \describe{
#' \item{strata_def}{Name for this strata scheme ('specht')}
#' \item{index}{Integer index, ordered by height level}
#' \item{label}{Stratum name (Ground, LowShrub, TallShrub, LowTree, MediumTree, TallTree)}
#' \item{lowerht}{Lower height of stratum (metres; -999.0 for ground layer)}
#' \item{upperht}{Upper height of stratum (metres; 999.0 for top layer)}
#' }
"StrataSpecht"


#' Strata definitions used by CERMB.
#'
#' Strata used within the Centre for Environmental Risk Management of Bushfires,
#' University of Wollongong, to rasterize LiDAR point clouds. The lowest two
#' strata are 0-0.3m and 0.3-0.5m. Subsequent strata are 0.5m intervals up to
#' 10m; then 1m intervals up to 30m; then 5m intervals up to 80m; then a
#' catch-all layer for >80m.
#'
#' @format A data frame with 6 rows and 5 columns:
#' \describe{
#' \item{strata_def}{Name for this strata scheme ('cermb')}
#' \item{index}{Integer index, ordered by height level}
#' \item{label}{Stratum label (Ground, 'upto0.5m', 'upto1m' ... 'over80m')}
#' \item{lowerht}{Lower height of stratum (metres)}
#' \item{upperht}{Upper height of stratum (metres)}
#' }
"StrataCERMB"


#' Former strata definitions used by CERMB.
#'
#' Strata formerly used within the Centre for Environmental Risk Management of
#' Bushfires, University of Wollongong, to rasterize LiDAR point clouds. These
#' have now been replaced by \code{\link{StrataCERMB}}. The lowest two strata
#' are 0-0.3m and 0.3-0.5m. Subsequent strata are 0.5m intervals up to 30m; then
#' a catch-all layer for >30m.
#'
#' @format A data frame with 6 rows and 5 columns:
#' \describe{
#' \item{strata_def}{Name for this strata scheme ('cermb_old')}
#' \item{index}{Integer index, ordered by height level}
#' \item{label}{Stratum label (Ground, 'upto0.5m', 'upto1m' ... 'over30m')}
#' \item{lowerht}{Lower height of stratum (metres)}
#' \item{upperht}{Upper height of stratum (metres)}
#' }
"StrataCERMB_OLD"
