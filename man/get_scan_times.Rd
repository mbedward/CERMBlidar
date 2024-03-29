% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tile_functions.R
\name{get_scan_times}
\alias{get_scan_times}
\title{Get scan times for a LAS tile}
\usage{
get_scan_times(las, by)
}
\arguments{
\item{las}{A LAS object, e.g. imported using \code{prepare_tile}.}

\item{by}{One of 'all' (default) or 'flightline'. Case-insensitive and may
be abbreviated.}
}
\value{
A data frame of start and end times.
}
\description{
The data table for an imported LAS tile includes a \code{gpstime} column
which gives, for each point, scan time expressed as \code{S - 1e9} where
\code{S} is the number of seconds since GPS epoch time: 1980-01-06 00:00:00
(GMT / UTC). This function converts the GPS time values to \code{POSIXct}
date-times and finds the start and end values, either for the tile as a whole
or for individual flightlines.
}
\seealso{
\code{\link{remove_zero_gpstimes}}
}
