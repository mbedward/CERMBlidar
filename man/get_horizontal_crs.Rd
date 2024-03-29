% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tile_functions.R
\name{get_horizontal_crs}
\alias{get_horizontal_crs}
\title{Get the horizontal coordinate reference system for a spatial object}
\usage{
get_horizontal_crs(x)
}
\arguments{
\item{x}{An object with an associated coordinate reference system that can
be retrieved using function \code{sf::st_crs()}. Usually this will be a
LAS object, e.g. imported using \code{prepare_tile}, but any other spatial
object class is acceptable.}
}
\value{
An \code{sf} package \code{crs} object.
}
\description{
Recent LiDAR data (e.g. from 2020 in Australia) can have a compound
coordinate reference system defined that consists of a horizontal and a
vertical component. Usually, we just want to use the horizontal component for
derived raster layers. This function extracts the horizontal component and
returns it as a \code{crs} object as used by the \code{sf} package. If the
input object has a simple (i.e. non-compound) coordinate reference system,
it will be returned unchanged.
}
