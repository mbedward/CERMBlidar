% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/flightline_functions.R
\name{get_flightline_info}
\alias{get_flightline_info}
\title{Get summary information for flight line extents}
\usage{
get_flightline_info(las, angular.tol = 25, min.ratio = 1.2, min.points = 1000)
}
\arguments{
\item{las}{A LAS object, e.g. imported using \code{prepare_tile}.}

\item{angular.tol}{The angular tolerance in degrees to apply when determining
flight line orientation. It refers to the angle between the longest side of
bounding rectangle and the horizontal (X coordinate axis). The default
value is 25 degrees and the valid range of values is
\code{0 < angular.tol <= 30}.}

\item{min.ratio}{The minimum ratio (default = 1.2) of the longest to the
shortest side of a flight line bounding rectangle for orientation to be
either 'NS' or 'EW'. Must be a value greater than 1.1. For flight lines
with more equilateral bounding rectangles, orientation is determined by
examining point GPS times (see Details).}

\item{min.points}{The minimum number of points in a flight line for it to be
included. Set to zero to include all flightlines. The default value (1000)
is intended to exclude flight lines that only appear at the margins of the
tile. Normally, such sparse flight lines will be removed when importing the
tile. In the case that all flight lines are included, bounding rectangles
and orientations will not be defined for those with very few points.}
}
\value{
A spatial data frame (class \code{sf}) with columns: flightlineID,
  xlen, ylen, orientation and geometry (minimum bounding rectangle of flight
  line).
}
\description{
This function determines the minimum bounding rectangle for points in each
flight line and, based on this, the orientation of the flight line: one of
'NS' (north - south), 'EW' (east - west), or XX (indeterminate). Note: this
function pre-dates the similar function \code{get_flightline_polygons} and
the two might be merged in a future version of the package.
}
\details{
To determine orientation, the function first checks the angle between the
longest side of the rectangle and the horizontal (X coordinate axis). If the
longest side is close to horizontal (where 'close' means plus or minus a
tolerance value specified with the \code{angular.tol} argument which has a
default value of 30 degrees) the orientation is provisionally labelled as
'EW', whereas if the longest side is close to vertical (Y coordinate axis) it
is provisionally labelled as 'NS'. Flight lines that are not sufficiently
close to horizontal or vertical are labelled as 'XX' and not considered
further.

For flight lines labelled 'NS' or 'EW', if the ratio of the longest side of
the bounding rectangle to the shortest side is greater than that specified by
the \code{ratio.side} argument (default 1.2), the initial orientation is
taken as confirmed. For flight lines with more equilateral bounding
rectangles, the point GPS times are checked along the X and Y dimensions to
determine the final orientation.
}
\seealso{
\code{\link{check_flightline_orientation}}
}
