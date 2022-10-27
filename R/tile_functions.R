#' Get names of mandatory input fields for LAS tiles
#'
#' Returns a named vector of character strings, where names are standard, single
#' character field abbreviations (as understood by \code{prepare_tile} and
#' \code{lidR::readLAS}) and values are field names. These are the fields that
#' must be included when importing a tile with \code{prepare_tile}.
#'
#' @return A named vector of lower case field names where names are standard,
#'   single character field abbreviations.
#'
#' @seealso \code{prepare_tile}
#'
#' @export
#'
mandatory_input_fields <- function() {
  c('x' = "x",
    'y' = "y",
    'z' = "z",
    't' = "gpstime",
    'c' = "classification")
}


#' Get bounds from the header of a LAS file or object
#'
#' This function reads the bounding coordinates and projection information from
#' the header of either a LAS file on disk (either compressed or uncompressed)
#' or LAS object imported via \code{prepare_tile}. The bounds are returned as
#' either: a named vector, an EWKT (Extended Well Known Text) string specifier,
#' an \code{sf} package bounding box, an \code{sf} package geometry list object
#' or a \code{terra} package \code{SpatExtent} object.
#'
#' @param x Either a character string specifying the path to a file, or a
#'   \code{LAS} object (e.g. imported using function \code{prepare_tile}).
#'   If a path, the file can be an uncompressed tile (\code{.las}), a
#'   compressed tile (\code{.laz}) or a zip file containing a \code{.las}
#'   file.
#'
#' @param type Either \code{'vec'} (default) to return a named vector of min and
#'   max coordinates and EPSG code; \code{'wkt'} to return an EWKT text string,
#'   \code{'sfc'} to return an \code{sf} package geometry list with a polygon,
#'   \code{'bbox'} to return an \code{sf} package bounding box object, or
#'   \code{'terra'} to return a \code{'terra'} package \code{'SpatExtent'}
#'   object.
#'
#' @param unzip.dir The directory in which to uncompress a compressed LAS file
#'   (identified by a '.zip' extension). If \code{NULL} (default) a temporary
#'   directory will be used. After processing, the uncompressed file will be
#'   deleted. Ignored if \code{x} is a LAS object or the path to a LAS or LAZ
#'   file.
#'
#' @return The bounding polygon in the format specified by the \code{type}
#'   argument.
#'
#' @examples
#' \dontrun{
#' # Read bounds from a zipped LAS file and return as a vector
#' bounds <- get_las_bounds("c:/somewhere/myfile.zip")
#'
#' # Get bounds as an 'sf' bounding box
#' bounds <- get_las_bounds("c:/somewhere/myfile.zip", type = "sf")
#' print(bounds)
#' plot(bounds)
#'
#' # Same thing for an imported LAS object
#' las <- prepare_tile("c:/somewhere/myfile.zip")
#' bounds <- get_las_bounds(las, type = "sf")
#' }
#'
#' @export
#'
get_las_bounds <- function(x, type = c("vec", "wkt", "sfc", "bbox", "terra"), unzip.dir = NULL) {
  type <- match.arg(type)

  if (inherits(x, "LAS")) {
    # Check for any missing coordinates
    if (anyNA(las@data$X) || anyNA(las@data$Y)) {
      stop("LAS tile has missing values for X and or Y. Fix this before proceeding.")
    }

    lims <-  c(xmin = min(las@data$X),
               xmax = max(las@data$X),
               ymin = min(las@data$Y),
               ymax = max(las@data$Y))

    xcrs <- get_horizontal_crs(x)

  } else if (is.character(x)) {
    if (length(x) > 1) {
      warning("Presently only one path element is supported. Ignoring extra elements.")
      x <- x[1]
    }

    if (!file.exists(x)) stop("File not found: ", x)

    unz.files <- character(0)
    if (.is_zipped(x)) {
      if (is.null(unzip.dir)) unzip.dir <- tempdir(check = TRUE)

      # allow for the possibility that there are other files in the zip file as well
      # as the LAS file
      unz.files <- utils::unzip(x, overwrite = TRUE, exdir = unzip.dir)

      las.file <- stringr::str_subset(unz.files, "\\.(las|LAS)$")
      n <- length(las.file)

      if (n == 0) {
        stop("zip file does not contain a file with extension las or LAS")
      }
      else if (n > 1) {
        stop("zip file contains multiple LAS files")
      }
    } else {
      las.file <- x
    }

    # Read header as list object
    hdr <- rlas::read.lasheader(las.file)
    if (length(unz.files) > 0) unlink(unz.files)

    lims <- c(xmin = hdr[["Min X"]],
              xmax = hdr[["Max X"]],
              ymin = hdr[["Min Y"]],
              ymax = hdr[["Max Y"]])

    # Convert header to LASheader object to get coord ref system
    hdr <- lidR::LASheader(hdr)

    suppressMessages(
      suppressWarnings({
        epsg <- lidR::epsg(hdr)

        if (epsg > 0) {
          xcrs <- sf::st_crs(epsg)
        } else {
          xcrs <- sf::st_crs(lidR::wkt(hdr))
        }
      })
    )
  } else {
    stop("Argument x should be a LAS object or the path to a LAS file")
  }


  if (type == "vec") {
    lims

  } else if (type == "terra") {
    terra::ext(lims)

  }  else {
    bb <- sf::st_bbox(lims, crs = xcrs)

    if (type == "bbox") {
      bb
    } else {
      glist <- sf::st_as_sfc(bb)
      if (type == "sfc") {
        glist
      } else {
        sf::st_as_text(glist, EWKT = TRUE)
      }
    }
  }
}


#' Get the horizontal coordinate reference system for a spatial object
#'
#' Recent LiDAR data (e.g. from 2020 in Australia) can have a compound
#' coordinate reference system defined that consists of a horizontal and a
#' vertical component. Usually, we just want to use the horizontal component for
#' derived raster layers. This function extracts the horizontal component and
#' returns it as a \code{crs} object as used by the \code{sf} package. If the
#' input object has a simple (i.e. non-compound) coordinate reference system,
#' it will be returned unchanged.
#'
#' @param x An object with an associated coordinate reference system that can
#'   be retrieved using function \code{sf::st_crs()}. Usually this will be a
#'   LAS object, e.g. imported using \code{prepare_tile}, but any other spatial
#'   object class is acceptable.
#'
#' @return An \code{sf} package \code{crs} object.
#'
#' @export
#'
get_horizontal_crs <- function(x) {
  xcrs <- sf::st_crs(x)
  if (is.na(xcrs)) return(xcrs)

  wkt <- sf::st_as_text(xcrs)

  if (!grepl("COMPD_CS", wkt)) {
    # Should just be a horizontal CRS - simply return it
    xcrs
  } else {
    # Extract the horizontal component
    i <- regexpr("PROJCS\\[", wkt)
    x <- substring(wkt, i)

    xc <- strsplit(x, "")[[1]]
    level <- 1
    k <- match("[", xc)
    while (level > 0) {
      k <- k + 1
      if (xc[k] == '[') {
        level <- level + 1
      } else if (xc[k] == ']') {
        level <- level - 1
      }
    }

    wkt <- substring(x, 1, k)

    sf::st_crs(wkt)
  }
}


#' Import and prepare a LAS tile for further processing
#'
#' This function imports data from a LAS file and prepares it for further
#' processing. Point heights are normalized to relative height above ground,
#' either by fitting a surface model to ground points (plus water points by
#' default) or by taking ground elevations from a provided raster DEM. Height
#' normalization can be disabled if desired. Flight lines are identified based
#' on GPS times for points and, optionally, any flight lines with less than a
#' threshold number of points are discarded. The flightlineID column is added
#' to the LAS header so that it will persist if the object is written to file.
#'
#' @param path Path to the LAS or LAZ format file to process. If the file
#'   extension is '.zip' it is assumed to be a compressed LAS file that will be
#'   unzipped before processing (see the \code{unzip.dir} parameter below). A
#'   compressed file should contain only one LAS file (identified by having a
#'   'las' or 'LAS' file extension) although it can also contain other files
#'   (e.g. HTML or XML documents).
#'
#' @param normalize.heights Whether, and how, to normalize point heights to
#'   ground level. Can be one of the following:
#'   \describe{
#'     \item{A logical value}{If \code{TRUE}, relative ground level is estimated from
#'       an elevation surface fitted by triangulation to ground points in the
#'       LAS file. This is equivalent to calling the
#'       \code{\link[lidR]{normalize}} function directly with the argument
#'       \code{algorithm = tin()}. If \code{FALSE}, point heights will not be
#'       normalized.}
#'     \item{An algorithm name as a character string}{Point heights will be
#'       normalized using the specified algorithm. Must be one of:
#'       \code{'tin', 'knnidw', 'kriging'} which correspond to the algorithm
#'       functions provided by the \code{lidR} package.}
#'     \item{A raster layer}{If a raster layer (a terra package \code{SpatRast}
#'       object or a raster package \code{RasterLayer} object) is provided, the
#'       cell values will be used as ground elevation to normalize point heights.
#'       The layer should have the same or greater extent as the LAS file.}
#'     \item{A raster filename}{Any character string that does not match one of
#'       the supported algorithm names will be treated as the path to a raster
#'       file. Supported formats are GeoTIFF ('.tif') and ESRI ASCII ('.asc').
#'       If the file extension is '.zip' it is assumed to be a compressed
#'       file that will be unzipped before processing (see the \code{unzip.dir}
#'       parameter below). A compressed file should contain only one raster file
#'       (identified by having a '.tif' or '.asc' file extension) although it
#'       can also contain other files (e.g. HTML or XML documents).}
#'     \item{NULL}{Same as \code{FALSE}, ie. point heights will not be
#'       normalized.}
#'   }
#'   \strong{The default value is \code{FALSE} to leave point heights
#'   unchanged.} If point heights are normalized, the original values are copied
#'   to a new data table column \code{'Zref'}, but this is not added to the
#'   LAS header by default so if the object is written to file the original
#'   heights will be lost.
#'
#' @param treat.as.ground Point classes to treat as ground points when
#'   normalizing point heights by interpolation using one of the lidR package
#'   algorithms (\code{'tin', 'knnidw', 'kriging'}). If point heights are being
#'   normalized from a raster DEM, any classes specified by this argument will
#'   have their normalized heights set to zero. The default value \code{c(2,9)}
#'   specifies ground and water classes.
#'
#' @param drop.negative If \code{TRUE} AND point heights are being normalized,
#'   any points whose heights are more than \code{drop.negative.threshold} below
#'   ground level will be discarded.
#'
#' @param drop.negative.threshold If \code{drop.negative} is \code{TRUE} AND
#'   point heights are being normalized, any heights more than this threshold
#'   value below ground level will be discarded. Default value is zero.
#'
#' @param drop.missing.coords If \code{TRUE} (default), discard any points that
#'   have missing values for the X, Y or Z coordinates. If \code{FALSE}, retain
#'   the point but issue a warning message.
#'
#' @param fields Either \code{NULL} (default) to include all data fields, or a
#'   character string containing single-letter abbreviations for selected
#'   fields. See \code{\link[lidR]{readLAS}} for details of the available,
#'   single-letter field abbreviations.
#'
#' @param classes Point classes to include or exclude. The default
#'   (\code{NULL}) is to include classes 2 (ground), 3-5 (vegetation), 6
#'   (buildings) and 9 (water). Specify a subset of classes as a vector of
#'   integers, e.g. \code{classes = 2:6} would include ground (2), vegetation
#'   (3, 4, 5) and building (6) points. Negative values can be used to exclude
#'   selected classes, e.g. \code{classes = -6} would include all classes except
#'   those classified as building points. Note that overlap points (class 12)
#'   are always excluded unless an explicit integer \code{classes} vector is
#'   provided that includes the value 12.
#'
#' @param min.points The minimum number of points in a flight line for it to be
#'   retained in the imported tile. The default value (1000) is intended to
#'   exclude flight lines that only appear at the margins of the tile.
#'
#' @param flight.gap The minimum time gap (seconds) to use when assigning points
#'   to flight lines.
#'
#' @param unzip.dir The directory in which to uncompress a compressed LAS file
#'   (identified by a '.zip' extension). If \code{NULL} (default) a temporary
#'   directory will be used. After processing, the uncompressed file will be
#'   deleted.
#'
#'
#' @return A \code{LAS} object.
#'
#' @export
#'
prepare_tile <- function(path,
                         normalize.heights = FALSE,
                         treat.as.ground = c(2,9),
                         drop.negative = TRUE,
                         drop.negative.threshold = 0,
                         drop.missing.coords = TRUE,
                         fields = NULL,
                         classes = NULL,
                         min.points = 1000,
                         flight.gap = 60,
                         unzip.dir = NULL) {

  # Flag that can be set below to indicate some change has been made to the
  # LAS data (e.g. dropping points) that requires the header to be updated.
  NEED_HEADER_UPDATE <- FALSE

  if (!is.character(path)) {
    stop("Argument path should be a character string for path and filename.")
  }

  if (!file.exists(path)) {
    stop("File not found: ", path)
  }

  if (length(path) > 1) {
    warning("Presently only one path element is supported. Ignoring extra elements.")
    path <- path[1]
  }

  if (is.null(fields)) fields <- "*"
  else fields <- .as_scalar(fields)

  treat.as.ground <- stats::na.omit(treat.as.ground)
  if (length(treat.as.ground) == 0) treat.as.ground <- c(2,9)

  # Vector to hold paths to (possibly) unzipped las and dem files
  unz.files <- character(0)

  zipped.las <- .is_zipped(path)
  if (zipped.las) {
    if (is.null(unzip.dir)) unzip.dir <- tempdir(check = TRUE)

    # allow for the possibility that there are other files in the zip file as well
    # as the LAS file
    ff <- utils::unzip(path, overwrite = TRUE, exdir = unzip.dir)

    las.file <- stringr::str_subset(ff, "\\.(las|LAS)$")
    n <- length(las.file)

    if (n == 0) {
      warning("zip file does not contain a file with extension las or LAS")
      return(NULL)
    }
    else if (n > 1) {
      warning("zip file contains multiple LAS files")
      return(NULL)
    }

    unz.files <- c(unz.files, ff)

  } else {
    las.file <- path
  }


  if (fields != "*") {
    # Check all mandatory fields are included and add any that were missed
    fields <- tolower(fields)
    must.haves <- names(mandatory_input_fields())
    found <- stringr::str_detect(fields, must.haves)
    if (any(!found)) {
      must.haves <- paste0(must.haves[!found], collapse="")
      fields <- paste0(fields, must.haves, collapse = "")
    }
  }

  if (is.null(classes)) {
    # default is ground (2), veg (3-5), buildings (6) and water (9)
    filtertxt <- "-keep_class 2 3 4 5 6 9"
  } else {
    keeps <- classes[ classes > 0 ]
    drops <- abs(classes[ classes < 0 ])

    # check for inconsistencies
    x <- base::intersect(keeps, drops)
    if (length(x) > 0) stop("Cannot keep and drop the same class(es): ", x)

    # Special treatment for class 12 (overlap points)
    if (!(12 %in% c(drops, keeps))) drops <- c(drops, 12)

    if (length(keeps) > 0)
      keepstxt <- paste("-keep_class", paste(keeps, collapse = " "))
    else
      keepstxt <- ""

    if (length(drops) > 0)
      dropstxt <- paste("-drop_class", paste(drops, collapse = " "))
    else
      dropstxt <- ""

    filtertxt <- paste(keepstxt, dropstxt)
  }

  las <- lidR::readLAS(las.file, select = fields, filter = filtertxt)

  # Check that the tile has some points (on rare occasions there are empty tiles)
  if (is_empty_tile(las)) {
    warning("Point cloud is empty in tile: ", path)
    return(las)
  }

  # Check for points with missing coordinates
  ii <- apply(las@data[, c("X", "Y", "Z")], MARGIN = 1, anyNA)
  if (any(ii)) {
    n <- sum(ii)
    if (drop.missing.coords) {
      message("Dropping ", n, " point(s) with missing coordinates")
      las@data <- las@data[!ii, ]
      NEED_HEADER_UPDATE <- TRUE

    } else {
      warning(n, " point(s) have missing coordinates but are not being dropped", immediate. = TRUE)
    }
  }

  # Normalize point heights relative to ground level if required
  # Check whether, and how, to normalize point heights
  do.normalize <-
    if (is.null(normalize.heights)) FALSE
    else if (is.logical(normalize.heights)) normalize.heights
    else TRUE

  if (do.normalize) {
    if (is.logical(normalize.heights)) { # must be TRUE
      if (normalize.heights) las <- lidR::normalize_height(las, algorithm = lidR::tin())

    } else if (is.character(normalize.heights)) {
      # Check if string is an algorithm name
      algorithms <- c("tin", "knnidw", "kriging")
      ii <- base::match(tolower(normalize.heights), algorithms)
      if (!is.na(ii)) {
        fn <- base::switch(algorithms[ii],
                           tin = lidR::tin,
                           knnidw = lidR::knnidw,
                           kriging = lidR::kriging)

        # Note the na.rm argument to avoid the process aborting if one or more points
        # cannot be normalized
        las <- lidR::normalize_height(las, algorithm = fn(), use_class = treat.as.ground, na.rm = TRUE)

      } else {
        # String should be a raster file path
        if (!file.exists(normalize.heights)) {
          stop("Cannot find raster file ", normalize.heights)
        }
        zipped.dem <- .is_zipped(normalize.heights)
        if (zipped.dem) {
          if (is.null(unzip.dir)) unzip.dir <- tempdir(check = TRUE)

          # allow for the possibility that there are other files in the zip file as well
          # as the raster file
          ff <- utils::unzip(normalize.heights, overwrite = TRUE, exdir = unzip.dir)

          dem.file <- stringr::str_subset(ff, "\\.(asc|ASC|tif|TIF)$")
          n <- length(dem.file)

          if (n == 0) {
            warning("DEM zip file does not contain a file with extension asc or tif")
            return(NULL)
          }
          else if (n > 1) {
            warning("DEM zip file contains multiple raster files")
            return(NULL)
          }

          unz.files <- unique(c(unz.files, ff))

        } else {
          dem.file <- normalize.heights
        }

        r <- terra::rast(dem.file)
        las <- .normalize_heights_dem(las, rdem = r)
      }

    } else if (.is_raster_object(normalize.heights) || .is_terra_object(normalize.heights)) {
      las <- .normalize_heights_dem(las, rdem = normalize.heights)

    } else {
      stop("Argument normalize.heights should be a logical value, ",
           "a character string, a SpatRast (terra package) or a ",
           "RasterLayer (raster package)")
    }

    # discard any points more than the threshold value below ground level
    if (drop.negative) {
      # Treat threshold value as height below zero
      ii <- las@data$Z < -abs(drop.negative.threshold)
      las@data <- las@data[!ii, ]
      NEED_HEADER_UPDATE <- TRUE
    }

    # Set the normalized height of all ground points to zero
    ii <- las@data$Classification %in% treat.as.ground
    las@data$Z[ii] <- 0
    NEED_HEADER_UPDATE <- TRUE
  }

  # Add flight line indices based on GPS times for points
  las <- lidR::retrieve_flightlines(las, dt = flight.gap)

  if (min.points > 0) las <- filter_flightlines(las, min.points)

  # Add the flightlineID column to the header
  las <- lidR::add_lasattribute(las, name = "flightlineID", desc = "Flight line ID")

  if (length(unz.files) > 0) unlink(unz.files)

  if (NEED_HEADER_UPDATE) {
    update_tile_header(las)
  } else {
    las
  }
}


.normalize_heights_dem <- function(las, rdem) {
  suppressWarnings(
    las <- mask_tile(las, rdem)
  )

  if (is.null(las)) {
    stop("No points lie within mask raster data cells")
  }

  lidR::normalize_height(las, algorithm = rdem)
}


#' Transform point cloud coordinates to a new projection
#'
#' This function takes the X,Y coordinates of all points and transforms them
#' from their existing projection to a user-specified one. The transformation is
#' done using \code{sf} package methods. If the point cloud is very large this
#' may take a while.
#'
#' @param las A LAS object, e.g. imported using \code{prepare_tile}.
#'
#' @param epsg.new Integer EPSG code for the projection into which the point
#'   cloud will be transformed.
#'
#' @param epsg.old This argument allows you to optionally supply or override the
#'   EPSG code for the existing projection. For a LAS tile with a properly
#'   assigned projection this defaults to the existing EPSG code.
#'
#' @return A new \code{LAS} object with reprojected point coordinates and
#'   an updated header.
#'
#' @export
#'
reproject_tile <- function(las, epsg.new, epsg.old = lidR::epsg(las)) {
  if (is.na(epsg.old) || is.null(epsg.old) || epsg.old <= 0) {
    stop("Existing projection is not defined. Specify a value for epsg.old")
  }

  if (epsg.new != epsg.old) {
    dat <- as.data.frame( las@data[, c("X", "Y")] )
    dat <- sf::st_as_sf(dat, coords = 1:2, crs = lidR::epsg(las))
    dat <- sf::st_transform(dat, epsg.new)
    dat <- sf::st_coordinates(dat)

    las@data$X <- dat[, "X"]
    las@data$Y <- dat[, "Y"]
    lidR::epsg(las) <- epsg.new

    las <- update_tile_header(las)
  }

  las
}


#' Apply a raster mask to a LAS tile
#'
#' This function uses a raster layer as a mask to subset a LAS point cloud.
#' Those points that fall within raster cells with non-missing values are
#' retained. Any points outside the bounds of the raster, or falling in cells
#' with \code{NA} values are dropped. If no points fall within data cells of
#' the raster, a warning message is issued and \code{NULL} is returned.
#'
#' @note This function \strong{does not check} that the coordinate reference
#'   systems of the LAS and raster objects are the same.
#'
#' @param las A LAS object, e.g. imported using \code{prepare_tile}.
#'
#' @param r A \code{RasterLayer} object to use as the mask.
#'
#' @return A new \code{LAS} object containing the masked points and updated
#'   header information.
#'
#' @export
#'
mask_tile <- function(las, r) {
  xy <- cbind(X = las@data$X, Y = las@data$Y)
  r <- .as_spat_raster(r)
  keep <- !is.na(terra::extract(r, xy))

  if (!any(keep)) {
    warning("No points fall within data cells of the raster mask")
    NULL
  } else {
    las@data <- las@data[keep, , drop = FALSE]
    update_tile_header(las)
  }
}


#' Filter flight lines based on the number of points in each
#'
#' This function takes an input LAS tile and returns a copy from which any
#' flight lines with less than a specified minimum number of points have been
#' removed. It is used by \code{prepare_tile} but can also be called
#' directly. If points are dropped, the header of the returned LAS object is
#' updated.
#'
#' @param las A LAS object, e.g. imported using \code{prepare_tile}.
#'
#' @param min.points The minimum number of points in a flight line for it to be
#'   retained. The default value (1e3) is intended to exclude flight lines that
#'   only appear at the margins of the tile.
#'
#' @return A (possibly) modified copy of the input tile with any sparse flight
#'   lines removed.
#'
#' @export
#'
filter_flightlines <- function(las, min.points = 1e3) {
  counts <- c( table(las@data$flightlineID) )
  ikeep <- counts >= min.points

  if (!any(ikeep)) {
    stop("No flightlines have the required minimum number of points")
  }

  if (all(ikeep)) {
    # nothing to do
  } else {
    # Identify points to keep and reduce las data to just
    # those records
    ids.keep <- as.integer( names(counts[ikeep]) )
    rec.keep <- las@data$flightlineID %in% ids.keep
    las@data <- las@data[rec.keep, ]

    # Update las header
    las <- update_tile_header(las)
  }

  las
}


#' Calculate average point density for a tile
#'
#' This function counts the points in the tile and divides by the total area
#' (based on the applicable distance units, e.g. square metres). For a tile with
#' average point density \code{p}, an estimate of the average distance between
#' neighbouring points, disregarding edge effects, is \code{0.5 / sqrt(p)}.
#'
#' @param las A LAS object, e.g. imported using \code{prepare_tile}.
#'
#' @return The average point density as a single value.
#'
#' @export
#'
get_point_density <- function(las) {
  np <- nrow(las@data)
  if (np == 0) stop("The tile is empty")

  dx <- diff(range(las@data$X))
  dy <- diff(range(las@data$Y))

  area <- dx * dy
  np / area
}


#' Find the minimum bounding rectangle that encloses specified points.
#'
#' When called with a single argument for the LAS tile, this function returns
#' the minimum bounding rectangle for the point cloud. The \code{classes} and
#' \code{flightlineIDs} arguments can be used to define a subset of points
#' to consider. If another kind of subset is required, the helper function
#' \code{link{.min_rectangle}} can be used directly.
#'
#' @param las A LAS object, e.g. imported using \code{prepare_tile}.
#'
#' @param classes Point classes to include: either a vector of integer class
#'   numbers or the string \code{'all'} (default).
#'
#' @param flightlines Flight lines to include: either a vector of integer
#'   class numbers or the string \code{'all'} (default).
#'
#' @return A data frame of corner vertices for the minimum bounding rectangle.
#'
#' Adapted from code by 'whuber' at: https://gis.stackexchange.com/a/22934/59514
#'
#' @export
#'
get_bounding_rectangle <- function(las, classes = "all", flightlines = "all") {
  if ("flightlineID" %in% colnames(las@data)) {

    lines.all <- sort(unique(las@data$flightlineID))

    if (tolower(flightlines) == "all") {
      flightlines <- lines.all
    } else if (length(flightlines) < 1) {
      stop("Argument flightlines should either be 'all' or a vector of one or more integer values")
    } else {
      ii <- flightlines %in% lines.all
      if (!all(ii)) {
        if (!any(ii)) stop("None of the specified flight lines appear in the tile")
        else {
          warning("Ignoring specified flight lines that are not in the tile")
          flightlines <- flightlines[ii]
        }
      }
    }
  }

  classes.all <- sort(unique(las@data$Classification))

  if (tolower(classes) == "all") {
    classes <- classes.all
  } else if (length(classes) < 1) {
    stop("Argument classes should either be 'all' or a vector of one or more integer values")
  } else {
    ii <- classes %in% classes.all
    if (!all(ii)) {
      if (!any(ii)) stop("None of the specified classes appear in the tile")
      else {
        warning("Ignoring specified classes that are not in the tile")
        classes <- classes[ii]
      }
    }
  }

  las@data %>%
    as.data.frame() %>%

    dplyr::filter(Classification %in% classes,
                  flightlineID %in% flightlines) %>%

    dplyr::select(X, Y) %>%

    .min_rectangle()
}


#' Find the minimum bounding rectangle that encloses a set of points
#'
#' This is a helper function called by \code{\link{get_bounding_rectangle}}.
#'
#' @param X,Y The X and Y arguments provide the point coordinates. Any
#'   reasonable way of defining the coordinates is acceptable. See the function
#'   \code{\link[grDevices]{xy.coords}} for details. If supplied separately,
#'   they must be of the same length.
#'
#' @return A data frame of corner vertices for the minimum bounding rectangle.
#'
#' Adapted from code by 'whuber' at: https://gis.stackexchange.com/a/22934/59514
#'
#' @export
#'
.min_rectangle <- function(X, Y = NULL) {
  xy <- grDevices::xy.coords(X, Y)

  # find vertices of the convex hull
  xy <- cbind(xy$x, xy$y)
  ii <- grDevices::chull(xy)
  m <- xy[ii, ]

  if (nrow(m) < 3) {
    warning("Minimum rectangle is undefined (e.g. single point or collinear points)")
    return(data.frame(X = NA_real_, Y = NA_real_))
  }

  # edge lengths
  dm <- cbind(
    dX = diff(m[,1]),
    dY = diff(m[,2])
  )

  # edge lengths
  lens <- apply(dm, 1, function(x) sqrt(x %*% x))

  # unit edge directions
  v <- diag(1/lens) %*% dm

  # normal directions to the edges
  w <- cbind(-v[,2], v[,1])

  # extremes along edges
  x <- apply(m %*% t(v), 2, range)
  y <- apply(m %*% t(w), 2, range)

  areas <- (y[1,]-y[2,])*(x[1,]-x[2,])
  k <- which.min(areas)

  rect <- cbind(x[c(1,2,2,1,1),k], y[c(1,1,2,2,1),k]) %*% rbind(v[k,], w[k,])
  colnames(rect) <- c("X", "Y")

  as.data.frame(rect)
}


#' Get scan times for a LAS tile
#'
#' The data table for an imported LAS tile includes a \code{gpstime} column
#' which gives, for each point, scan time expressed as \code{S - 1e9} where
#' \code{S} is the number of seconds since GPS epoch time: 1980-01-06 00:00:00
#' (GMT / UTC). This function converts the GPS time values to \code{POSIXct}
#' date-times and finds the start and end values, either for the tile as a whole
#' or for individual flightlines.
#'
#' @param las A LAS object, e.g. imported using \code{prepare_tile}.
#'
#' @param by One of 'all' (default) or 'flightline'. Case-insensitive and may
#'   be abbreviated.
#'
#' @return A data frame of start and end times.
#'
#' @seealso \code{\link{remove_zero_gpstimes}}
#'
#' @export
#'
get_scan_times <- function(las, by) {
  if (missing(by) || is.null(by) || by == "") {
    by <- "all"
  } else {
    by <- match.arg(tolower(by), c("all", "flightline"))
  }

  times <- .get_point_gpstimes(las)

  if (by == "all") {
    data.frame(time.start = min(times), time.end = max(times))

  } else {
    data.frame(flightlineID = las@data$flightlineID, times) %>%
      dplyr::group_by(flightlineID) %>%
      dplyr::summarize(time.start = min(times), time.end = max(times)) %>%
      dplyr::arrange(flightlineID)
  }
}


#' Check for and drop any points with zero GPS times
#'
#' Sometimes the GPS time for a point is zero rather than the expected large
#' integer value (seconds since epoch), which will lead to incorrect capture
#' times being reported for the LiDAR tile. This function checks for the
#' presence of such invalid GPS times and drops those points from the data if
#' they constitute no more than a specified proportion of the total points. The
#' function returns the (possibly modified) LAS object.
#'
#' Note: Any points with a missing value for GPS time will be retained.
#'
#' @param las A LAS object, e.g. imported using \code{prepare_tile}.
#'
#' @param max_prop The maximum proportion of points with zero GPS times. This
#'   should be set to a small value (default is 0.01). If the proportion of zero
#'   GPS times exceeds this value, the function will throw an error.
#'
#' @param quiet Set to \code{TRUE} to suppress messages. Default is
#'   \code{FALSE}.
#'
#' @return The (possibly modified) LAS object.
#'
#' @examples
#' \dontrun{
#' las <- remove_zero_gpstimes(las)
#' }
#'
#' @seealso \code{\link{get_scan_times}}
#'
#' @export
#'
remove_zero_gpstimes <- function(las, max_prop = 0.01, quiet = FALSE) {
  # Using `which` here will safely ignore any NA values in gpstime
  i_zero <- which(las@data$gpstime == 0)
  n_zero <- length(i_zero)

  if (n_zero > 0) {
    p <- n_zero / nrow(las@data)

    if (p <= max_prop) {
      if (!quiet) {
        msg <- glue::glue("Removing {nzero} points with zero GPS time")
        message(msg)
      }

      las@data <- las@data[-i_zero, ]
      las <- update_tile_header(las)

    } else {
      msg <- glue::glue("Proportion of points with zero GPS time ({round(p, 4)})
                         is greater than the maximum allowable threshold ({max_prop})")
      stop(msg)
    }
  } else {
    if (!quiet) {
      message("No points have zero GPS time")
    }
  }

  las
}


# Private helper function to get point GPStimes as POSIXct
# date-times
#
.get_point_gpstimes <- function(las) {
  T0 <- lubridate::ymd_hms("1980-01-06 00:00:00", tz = "GMT")
  times <- T0 + lubridate::seconds(1e9 + las@data$gpstime)

  times
}


#' Check whether a LAS tile object has been prepared
#'
#' When a LAS tile is imported using the \code{prepare_tile} function,
#' point heights are normalized relative to ground elevation and flight lines
#' are identified. This results in two extra columns (Zref and flightlineID)
#' being added to the LAS data table. This function checks that these two
#' columns are present. It also checks that the tile has point data, returning
#' FALSE if the point cloud is empty.
#'
#' @param las A LAS object, e.g. imported using \code{prepare_tile}.
#'
#' @return TRUE if the tile has been prepared or FALSE otherwise.
#'
#' @export
#'
is_prepared_tile <- function(las) {
  if (!inherits(las, "LAS")) stop("Object is not a LAS tile")

  if (is_empty_tile(las)) {
    FALSE
  } else {
    expected.cols <- c(mandatory_input_fields(), "zref", "flightlineid")

    all(expected.cols %in% tolower(colnames(las@data)))
  }
}



#' Check whether a LAS tile has an empty point cloud
#'
#' Checks whether a point cloud is empty, i.e. there are zero rows in the LAS
#' object's data table.
#'
#' @param las A LAS object, e.g. imported using \code{prepare_tile}.
#'
#' @return TRUE if the point cloud is empty.
#'
#' @export
#'
is_empty_tile <- function(las) {
  if (!inherits(las, "LAS")) stop("Object is not a LAS tile")
  nrow(las@data) == 0
}


#' Create raster layers summarizing point counts for a LAS tile
#'
#' Given a LAS tile and a raster cell size, this function creates a raster
#' object with a band for each point class being considered. Cell values in each
#' band give the point count over all heights for the given point class. A layer
#' can also be created for all three vegetation classes combined. The resulting
#' layers are returned as a \code{terra::SpatRast} object.
#'
#' @param las A LAS object, e.g. imported using \code{prepare_tile}.
#'
#' @param res Raster cell size. The default value is 10 which assumes that
#'   the LAS data is in a projected coordinate system with metres as units.
#'
#' @param classes An integer vector specifying the point classes to summarize.
#'   May include 2 (ground), 3 (low veg), 4 (mid veg), 5 (high veg), 35 (all veg
#'   - the three vegetation classes combined), 6 (buildings) and 9 (water). The
#'   default is to summarize all of the above.
#'
#' @return A \code{terra::SpatRast} object with a layer for each requested
#'   attribute.
#'
#' @examples
#' \dontrun{
#' # Get summary counts at 10m raster cell size for all supported
#' # point classes
#' rcounts <- get_summary_counts(las)
#'
#' # Get summary counts for ground and combined vegetation classes
#' # at 2m resolution
#' rcounts <- get_summary_counts(las, res = 2, classes = c("ground", "allveg"))
#' }
#'
#' @export
#'
get_summary_counts <- function(las, res = 10, classes = NULL) {

  ClassLookup <- data.frame(
    code = c(2, 3, 4, 5, 35, 6, 9),
    label = c('ground',
              'low veg', 'mid veg', 'high veg', 'all veg',
              'building', 'water')
  )

  if (is.null(classes) || length(classes) == 0) {
    classes <- ClassLookup$code

  } else if (is.numeric(classes)) {
    classes <- unique(classes)
    ok <- classes %in% ClassLookup$code
    if (!all(ok)) {
      stop("Unrecognized class(es): ", classes[!ok], "\n",
           "Valid options: ", paste(ClassLookup$code, collapse = ", "))
    }

  } else {
    stop("Argument classes should be NULL or a vector of one or more integers")
  }

  ext <- get_las_bounds(las, type = "terra")
  rbase <- terra::rast(ext, resolution = res)

  rr <- lapply(classes, function(cl) {
    if (cl == 35) cl <- 3:5
    ii <- las@data$Classification %in% cl

    if (!any(ii)) {
      warning("No points in class ", cl)
      r <- terra::init(rbase, 0)
    } else {
      p <- terra::vect( as.matrix(las@data[ii, c("X", "Y"), drop=FALSE]) )
      r <- terra::rasterize(p, rbase, fun = length, background = 0)
    }

    r
  })

  res <- terra::rast(rr)
  ii <- match(classes, ClassLookup$code)
  names(res) <- ClassLookup$label[ii]

  terra::crs(res) <- sf::st_as_text( get_horizontal_crs(las) )

  res
}


#' Get LiDAR point counts for defined strata.
#'
#' Give a LAS tile and a set of height ranges, this function creates a
#' multi-band raster where the cells in each band represent voxels for a given
#' stratum, and each cell value is the number of points within the voxel. This
#' can be used as the first step in the calculation of vegetation cover
#' estimates.
#'
#' @param las A LAS object, e.g. imported using \code{prepare_tile}.
#'
#' @param strata A data frame of strata definitions with (at least) three
#'   columns: label, lowerht, upperht. Strata must cover all heights and be
#'   non-overlapping. This will be checked by the function \code{check_strata}.
#'
#' @param res Raster cell size.
#'
#' @param classes Integer classification codes for points to include.
#'   Default is 2 (ground) and 3:5 (vegetation) and 9 (water).
#'
#' @return A \code{terra::SpatRast} object where each band corresponds to a
#'   stratum (including the ground layer) and cell values are voxel point
#'   counts.
#'
#' @seealso \code{\link{check_strata}}
#'
#' @export
#'
get_stratum_counts <- function(las, strata, res = 10, classes = c(2,3,4,5,9)) {
  if (!is_prepared_tile(las)) {
    message("No Zref column. Assuming column Z is normalized heights")
  }
  strata <- check_strata(strata)

  # Convert lower and upper strata heights to a set of cut-points
  # for use with the findInterval function.
  # This code assumes all upper heights are unique which should have
  # been checked by check_strata().
  #
  n <- nrow(strata)
  b <- unique(sort(c(strata$lowerht, strata$upperht)))

  stratabreaks = list(
    nstrata = n,
    breaks = b,
    strata.indices = match(b, strata$upperht)[-1]
  )

  xy <- cbind(X = lidR:::f_grid(las@data$X, res, 0),
              Y = lidR:::f_grid(las@data$Y, res, 0) )


  ext <- get_las_bounds(las, type = "terra")
  rbase <- terra::rast(ext, resolution = res)

  # points in desired classes
  include <- las@data$Classification %in% classes

  xy <- xy[include, , drop = FALSE]

  bins <- findInterval(las@data$Z[include],
                       stratabreaks$breaks,
                       left.open = TRUE)

  bins[bins == 0] <- NA

  xy <- cbind(xy, stratum = stratabreaks$strata.indices[bins])
  xy <- xy[ !is.na(bins), , drop = FALSE ]

  rcounts <- lapply(
    1:nrow(strata),
    function(i) {
      dat <- xy[ xy[,3] == i, , drop=FALSE]
      if (nrow(dat) == 0) {
        terra::init(rbase, 0)
      } else {
        p <- terra::vect( dat[, 1:2, drop = FALSE] )
        terra::rasterize(p, rbase, fun = length, background = 0)
      }
    })


  rcounts <- terra::rast(rcounts)
  names(rcounts) <- strata$label

  # Note: guard against compound CRS (horizontal + vertical) now used
  # by recent LiDAR files
  terra::crs(rcounts) <- sf::st_as_text( get_horizontal_crs(las) )

  rcounts
}


#' Extract building points from a LAS tile
#'
#' @param las A LAS object, e.g. imported using \code{prepare_tile}.
#'
#' @return An \code{sf} spatial data frame suitable for map display and export
#'   to a shapefile.
#'
#' @examples
#' \dontrun{
#' # Import tile and normalize point heights
#' las <- prepare_tile(path.to.tile)
#'
#' # Extract building points
#' buildings <- get_building_points(las)
#'
#' # Write to shapefile (requires 'sf' package)
#' st_write(buildings, "buildings.shp")
#' }
#'
#' @export
#'
get_building_points <- function(las) {
  if (!is_prepared_tile(las)) {
    message("No Zref column. Assuming column Z is normalized heights")
  }

  dat <- las@data %>%
    as.data.frame() %>%
    dplyr::filter(Classification == 6) %>%
    dplyr::select(X, Y, Z)

  sf::st_as_sf(dat, coords = c("X", "Y"),
               crs = get_horizontal_crs(las))
}


#' @title Calculate vegetation cover estimates for strata
#'
#' @description This function was formerly used to calculate maximum likelihood
#'   mean cover estimates for strata. Please now use function
#'   \code{\link{get_cover_quantiles}} instead with argument \code{mle=TRUE}.
#'
#' @param rcounts A raster object with a band of point count data for each
#'   vegetation stratum plus a band for the ground layer. Normally, this will be
#'   a \code{terra} package \code{SpatRaster} object generated by function
#'   \code{get_stratum_counts()}, but it may also be a \code{RasterStack} or
#'   \code{RasterBrick} (\code{raster} package) object.
#'
#' @param background Value (default \code{NA}) to assign as the cover value for
#'   a stratum when there are no LiDAR points for the stratum, any lower
#'   strata or ground.
#'
#' @return A \code{SpatRast} object with a band for each vegetation stratum.
#'
#' @templateVar fun get_stratum_cover
#' @template template-depr_fun
NULL

#' @templateVar old get_stratum_cover
#' @templateVar new get_cover_quantiles
#' @template template-depr_pkg
#'
#' @export
#'
get_stratum_cover <- function(rcounts, background = NA) {
  .Deprecated("get_cover_quantiles")

  get_cover_quantiles(rcounts, probs = NULL, mle = TRUE)
}


#' Derive a raster of maximum point height in all or selected classes
#'
#' This is useful for checking purposes. It can also be used to derive
#' an approximate vegetation height layer by setting the \code{classes}
#' argument to \code{c(3:5)}.
#'
#' @param las A LAS object, e.g. imported using \code{prepare_tile}.
#'
#' @param res Raster cell size. The default value is 10 which assumes that
#'   the LAS data is in a projected coordinate system with metres as units.
#'
#' @param classes Point classes to include when calculating maximum height
#'   for each cell. Default is c(2, 3:5, 9) for ground, vegetation and
#'   water point classes.
#'
#' @param point_height_range A two-element numeric vector with the minimum and maximum
#'   heights of points to include when determining the maximum height. The values can be
#'   in any order. The default is \code{c(0, Inf)}. This can be used, for example when
#'   determining approximate vegetation height, to filter out points with unusually high
#'   elevations (e.g. birds or airborne dust). See example below.
#'
#' @param nodata Numeric value for cells with no points. Default is \code{NA_real_}.
#'
#' @return A \code{RasterLayer} in which cell values are maximum point height.
#'
#' @examples
#' \dontrun{
#' # Determine approximate vegetation height for each voxel, setting
#' # the height range filter based on tallest New South Wales tree height
#' #
#' rmaxht <- get_max_height(las, classes = 3:5, point_height_range = c(0, 76))
#' }
#'
#' @export
#'
get_max_height <- function(las,
                           res = 10,
                           classes = c(2,3:5,9),
                           point_height_range = c(0, Inf),
                           nodata = NA_real_) {

  ok <- is.vector(point_height_range) &&
    is.numeric(point_height_range) &&
    (length(point_height_range) == 2)

  if (!ok) {
    stop("point_height_range should be a two-element numeric vector")
  }

  point_height_range <- sort(point_height_range)

  if (!is.na(nodata)) {
    if (is.null(nodata)) {
      nodata <- NA_real_
    } else if (!is.numeric(nodata)) {
      stop("Argument nodata should be NA or a numeric value")
    }
  }

  if (length(classes) == 0 || !is.numeric(classes)) {
    stop("Argument classes should be a vector of one or more integer values")
  }

  stopifnot(res > 0)

  las_sub <- lidR::filter_poi(las,
                              Classification %in% classes,
                              Z >= point_height_range[1],
                              Z <= point_height_range[2])

  #r <- lidR::grid_canopy(las_sub, res = res, algorithm = lidR::p2r(1))

  xy <- as.matrix(las_sub@data[, c("X", "Y")])
  z <- las_sub@data$Z

  ext <- get_las_bounds(las, type = "terra")
  rbase <- terra::rast(ext, resolution = res)

  # Still using raster function here because it is much faster than
  # creating a SpatVect point object and calling terra::rasterize
  r <- raster::rasterize(xy, raster::raster(rbase), field = z, fun = max, background = nodata)

  r <- terra::rast(r)
  terra::crs(r) <- sf::st_as_text( get_horizontal_crs(las) )

  r
}


#' Get point class frequencies.
#'
#' Gets the number of points in classes: ground, veg, building, water, and
#' other.
#'
#' @param las A LAS object.
#'
#' @return A named list of the count of points in each of the following
#' classes: ground, veg, building, water, other.
#'
#' @export
#'
get_class_frequencies <- function(las) {
  x <- table(las@data$Classification)
  x <- data.frame(code = names(x), n = as.numeric(x), stringsAsFactors = FALSE)

  x$class <- sapply(x$code, switch,
                    '2' = "ground",
                    '3' = "veg",
                    '4' = "veg",
                    '5' = "veg",
                    '6' = "building",
                    '9' = "water",
                    "other")

  x <- tapply(x$n, x$class, sum)

  all.classes <- c("ground", "veg", "building", "water", "other")
  to.add <- setdiff(all.classes, names(x))
  n <- length(to.add)
  if (n > 0) {
    x2 <- integer(n)
    names(x2) <- to.add
    x <- c(x, x2)
  }

  as.list(x)
}


#' Check strata lookup table and put it into a standard form
#'
#' This function is used by \code{\link{get_stratum_counts}}. It checks that:
#' \itemize{
#'   \item The object is a data frame with (at least) three columns labelled:
#'     label, lowerht, upperht. Column order and case of names does not matter.
#'   \item All stratum upper heights are greater than lower heights.
#'   \item No two strata overlap (gaps between strata are permitted).
#' }
#' An error results if any of these conditions are not met.
#'
#' @param strata A data frame of strata definitions.
#'
#' @return A data frame in standard form with (at least) columns: label,
#'   lowerht, upperht, and rows ordered by ascending stratum heights.
#'
#' @export
#'
check_strata <- function(strata) {
  if (!inherits(strata, "data.frame"))
    stop("Strata lookup table should be a data frame")

  names <- tolower( colnames(strata) )

  if (!all( c("label", "lowerht", "upperht") %in% names ))
    stop("Strata lookup table should have columns: label, lowerht, upperht")

  colnames(strata) <- names

  if (!all(strata$upperht > strata$lowerht))
    stop("All upper heights should be greater than lower heights ",
         "in strata lookup table")

  o <- order(strata$lowerht)
  strata <- strata[o, ]

  n <- nrow(strata)
  if (!all(strata$lowerht[-1] >= strata$upperht[-n]))
    stop("Strata lookup table has overlapping strata")

  strata
}


#' Create a transformation table for two sets of strata
#'
#' Creates a transformation table that can be used to aggregate a set of input
#' strata into a reduced number of broader output strata. The table can be used
#' with function \code{aggregate_stratum_counts} to combine the bands of a point
#' count raster. Usually the reduced strata will cover the same range of heights
#' as the input strata, but it is also possible to create a subset table, e.g.
#' only covering lower heights. The height boundaries of input and output
#' strata must align. For example, if the input strata included two levels for
#' 2-5m and 5-10m, these could be aggregated into a layer for 2-10m but not into
#' a layer for 2-8m.
#'
#' @param from_strata Data frame defining the input strata to be aggregated, in
#'   the format used by functions \code{check_strata} and
#'   \code{get_stratum_counts}. The number of strata should be greater than that
#'   of \code{to_strata}.
#'
#' @param to_strata Data frame defining the broader, transformed output strata,
#'   in the format used by functions \code{check_strata} and
#'   \code{get_stratum_counts}. The number of strata should be greater than that
#'   of \code{to_strata}.
#'
#' @param flag_lowerht Value(s) used for the lower height of the lower-most
#'   stratum, which is usually treated as having no lower bound. Default is
#'   \code{c(-Inf, -999.0)}
#'
#' @param flag_upperht Value(s) used for the upper height of the upper-most
#'   stratum, which is usually treated as having no upper bound. Default is
#'   \code{c(-Inf, -999.0)}
#'
#' @return A data frame defining the transformation.
#'
#' @seealso \code{\link{check_strata}} and \code{\link{get_stratum_counts}}
#'
#' @examples
#' xtrans <- get_strata_transform(CERMBlidar::StrataCERMB, CERMBlidar::StrataSpecht)
#'
#' # This transform is invalid
#'
#' @export
#'
get_strata_transform <- function(from_strata, to_strata,
                                 flag_lowerht = c(-Inf, -999.0),
                                 flag_upperht = c(Inf, 999.0)) {

  from_strata <- check_strata(from_strata)
  to_strata <- check_strata(to_strata)

  if (nrow(from_strata) <= nrow(to_strata)) {
    stop("from_strata should have more layers than to_strata")
  }

  # Check if transform will be valid
  ii <- !(to_strata$lowerht %in% flag_lowerht)
  ok <- all(to_strata$lowerht[ii] %in% from_strata$lowerht)
  if (!ok) stop("One or more lower heights in to_strata do not align with from_strata")

  ii <- !(to_strata$upperht %in% flag_upperht)
  ok <- all(to_strata$upperht[ii] %in% from_strata$upperht)
  if (!ok) stop("One or more upper heights in to_strata do not align with from_strata")

  # Create and return the transform table
  from_strata %>%
    dplyr::rename_with( ~ paste0("from_", .x) ) %>%

    dplyr::mutate(to_index = cut(from_upperht,
                          breaks = c(-Inf, to_strata$upperht),
                          labels = FALSE,
                          right = TRUE),

           to_label = to_strata$label[to_index] )
}


#' Update a LAS header based on current point data
#'
#' This function can be called to update the header information of a LAS object
#' after making direct changes to its points data (e.g. removing selected
#' points). It does a little more than the similar function \code{header_update}
#' in the \code{rlas} package which only updates the header elements for min and
#' max X, Y and Z, and the number of points. This function additionally updates
#' the X and Y offsets.
#'
#' @param las A LAS object, e.g. imported using \code{prepare_tile}.
#'
#' @return A LAS object with header information updated, if necessary, to
#'   accord with its current data.
#'
#' @export
#'
update_tile_header <- function(las) {
  hdr <- as.list(las@header)

  hdr[["Number of point records"]] <- nrow(las@data)

  if ("ReturnNumber" %in% colnames(las@data))
    hdr[["Number of points by return"]] <- as.numeric(table(las@data$ReturnNumber))

  xlims <- range(las@data$X)
  ylims <- range(las@data$Y)
  zlims <- range(las@data$Z)

  hdr[["Min X"]] <- xlims[1]
  hdr[["Min Y"]] <- ylims[1]
  hdr[["Min Z"]] <- zlims[1]

  hdr[["Max X"]] <- xlims[2]
  hdr[["Max Y"]] <- ylims[2]
  hdr[["Max Z"]] <- zlims[2]

  hdr[["X offset"]] <- round(xlims[1], -floor(log10(xlims[1])))
  hdr[["Y offset"]] <- round(ylims[1], -floor(log10(ylims[1])))

  las@header <- lidR::LASheader(hdr)
  las
}
