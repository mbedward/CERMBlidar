#' Import and prepare a LAS tile for further processing
#'
#' This function imports data from a LAS file and prepares it for further
#' processing. A surface model is first fitted to ground points (class 2) by
#' Delaunay triangulation and the elevation of all points is then adjusted to be
#' relative to ground level.
#'
#' @param path Path to the LAS tile to process. If the file extension is '.zip'
#'   it is assumed to be a compressed LAS file that will be unzipped before
#'   processing (see the \code{unzip.dir} parameter below). A compressed file
#'   should contain only one LAS file (identified by having a 'las' or 'LAS'
#'   file extension) although it can also contain other files (e.g. HTML or XML
#'   documents).
#'
#' @param normalize.heights If TRUE (default), normalize all point heights to
#'   ground level as estimated from an elevation surface fitted by Delaunay
#'   triangulation. The original point heights are kept as a new data table
#'   column 'Zref'.
#'
#' @param drop.negative If TRUE, any points (other than ground points) whose
#'   heights are below ground level (as estimated by Delaunay interpolation of
#'   ground point heights) are discarded. This argument only applies if
#'   \code{normalize.heights = TRUE}.
#'
#' @param fields A string containing single-letter abbreviations for the data
#'   fields to include. The default is to return all fields. See
#'   \code{\link[lidR]{readLAS}} for details of the available abbreviations.
#'
#' @param classes Point classifcations to include. Default is ground (2), vegetation
#'   (3, 4, 5) and buildings (6).
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
                         normalize.heights = TRUE,
                         drop.negative = TRUE,
                         fields = "*",
                         classes = c(2, 3, 4, 5, 6),
                         min.points = 1000,
                         flight.gap = 60,
                         unzip.dir = NULL) {

  if (length(path) > 1) {
    warning("Presently only one path element is supported. Ignoring extra elements.")
    path <- path[1]
  }

  zipped <- .is_zipped(path)
  if (zipped) {
    if (is.null(unzip.dir)) unzip.dir <- tempdir(check = TRUE)

    # allow for the possibility that there are other files in the zip file as well
    # as the LAS file
    unz.files <- utils::unzip(path, overwrite = TRUE, exdir = unzip.dir)

    las.file <- stringr::str_subset(unz.files, "\\.(las|LAS)$")
    n <- length(las.file)

    if (n == 0) stop("zip file does not contain a file with extension las or LAS")
    if (n > 1) stop("zip file contains multiple LAS files")

  } else {
    las.file <- path
  }

  filtertxt <- paste("-keep_class", paste(classes, collapse = " "))

  las <- lidR::readLAS(las.file, select = fields, filter = filtertxt)

  # Check that the tile has some points (rarely we encouter empty tiles)
  if (nrow(las@data) == 0) {
    stop("Point cloud is empty in tile read from ", path)
  }

  # Normalize point heights relative to ground level
  if (normalize.heights) {
    las <- lidR::lasnormalize(las, algorithm = lidR::tin())

    # remove negative heights
    if (drop.negative) {
      ii <- las@data$Classification == 2 | las@data$Z > 0
      las@data <- las@data[ ii, ]
    }
  }

  # Add flight line indices based on GPS times for points
  #lidR::lasflightline(dt = flight.gap)
  las <- .tag_flight_line(las, dt = flight.gap)

  if (min.points > 0) las <- filter_flightlines(las, min.points)

  if (zipped) {
    unlink(unz.files)
  }

  las
}


# This function is called by prepare_tile to identify flight lines based
# on point GPS times. The function lidR::lasflightline does this but it
# was broken in version 2.0.1 of the package, hence this work-around.
.tag_flight_line <- function(las, dt = 30) {
  x <- data.frame(t = las@data[, "gpstime"])
  x$i <- 1:nrow(x)
  x <- x[order(x$t),]

  b <- c(TRUE, diff(x$t) > dt)
  b[1] <- TRUE
  x$flightlineID <- cumsum(b)

  x <- x[order(x$i),]

  las@data$flightlineID <- x$flightlineID

  las
}


.is_zipped <- function(path) stringr::str_detect(path, "\\.zip\\s*$")


#' Filter flight lines based on the number of points in each
#'
#' This function takes an input LAS tile and returns a copy from which any
#' flight lines with less than a specified minimum number of points have been
#' removed. It is used by \code{\link{prepare_tile}} but can also be called
#' directly.
#'
#' @param las A LAS object, e.g. imported using \code{\link{prepare_tile}}.
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
  ii <- counts >= min.points

  if (!any(ii)) {
    stop("No flightlines have the required minimum number of points")
  }

  ids.keep <- as.integer( names(counts[ii]) )

  irecs <- las@data$flightlineID %in% ids.keep
  las@data <- las@data[irecs, ]

  las
}

#' Check that the orientations of all flightlines match
#'
#' This is a convenience function that calls \code{\link{get_flightline_dir}}
#' with default arguments and checks that all flight lines have the same
#' orientation, either 'EW' or 'NS'.
#'
#' @param las A LAS object, e.g. imported using \code{\link{prepare_tile}}.
#'
#' @param ... Additional arguments to pass to \code{\link{get_flightline_info}}
#'   including \code{angular.tol}.
#'
#' @return \code{TRUE} if flight lines are consistent; \code{FALSE} otherwise.
#'
#' @seealso \code{\link{get_flightline_info}}
#'
#' @export
#'
check_flightlines <- function(las, ...) {
  dat <- get_flightline_info(las, ...)
  o <- na.omit(dat$orientation)

  if (length(o) == 0) FALSE
  else (o[1] %in% c("EW", "NS")) && all(o == o[1])
}


#' Get summary information for flight line extents
#'
#' This function determines the minimum bounding rectangle for points in each
#' flight line and, based on this, the orientation of the flight line: one of
#' 'NS' (north - south), 'EW' (east - west), or XX (indeterminate).
#'
#' To determine orientation, the function first checks the angle between the
#' longest side of the rectangle and the horizontal (X coordinate axis). If the
#' longest side is close to horizontal (where 'close' means plus or minus a
#' tolerance value specified with the \code{angular.tol} argument which has a
#' default value of 30 degrees) the orientation is provisionally labelled as
#' 'EW', whereas if the longest side is close to vertical (Y coordinate axis) it
#' is provisionally labelled as 'NS'. Flight lines that are not sufficiently
#' close to horizontal or vertical are labelled as 'XX' and not considered
#' further.
#'
#' For flight lines labelled 'NS' or 'EW', if the ratio of the longest side of
#' the bounding rectangle to the shortest side is greater than that specified by
#' the \code{ratio.side} argument (default 1.2), the initial orientation is
#' taken as confirmed. For flight lines with more equilateral bounding
#' rectangles, the point GPS times are checked along the X and Y dimensions to
#' determine the final orientation.
#'
#'
#' @param las A LAS object, e.g. imported using \code{\link{prepare_tile}}.
#'
#' @param angular.tol The angular tolerance in degrees to apply when determining
#'   flight line orientation. It refers to the angle between the longest side of
#'   bounding rectangle and the horizontal (X coordinate axis). The default
#'   value is 25 degrees and the valid range of values is
#'   \code{0 < angular.tol <= 30}.
#'
#' @param min.ratio The minimum ratio (default = 1.2) of the longest to the
#'   shortest side of a flight line bounding rectangle for orientation to be
#'   either 'NS' or 'EW'. Must be a value greater than 1.1. For flight lines
#'   with more equilateral bounding rectangles, orientation is determined by
#'   examining point GPS times (see Details).
#'
#' @param min.points The minimum number of points in a flight line for it to be
#'   included. Set to zero to include all flightlines. The default value (1000)
#'   is intended to exclude flight lines that only appear at the margins of the
#'   tile. Normally, such sparse flight lines will be removed when importing the
#'   tile. In the case that all flight lines are included, bounding rectangles
#'   and orientations will not be defined for those with very few points.
#'
#' @return A spatial data frame (class \code{sf}) with columns: flightlineID,
#'   xlen, ylen, orientation and geometry (minimum bounding rectangle of flight
#'   line).
#'
#' @importFrom dplyr %>%
#'
#' @seealso \code{\link{check_flightlines}}
#'
#' @export
#'
get_flightline_info <- function(las,
                                angular.tol = 25,
                                min.ratio = 1.2,
                                min.points = 1000) {

  stopifnot(angular.tol > 0.0, angular.tol <= 30)
  stopifnot(min.ratio > 1.1)

  ids <- sort(unique(las@data$flightlineID))

  # Count points in each flight line.
  # The c() call is to convert the table object to a
  # simple vector.
  counts.all <- c( table(las@data$flightlineID) )
  total.count <- sum(counts.all)

  # Check for flightlines with an insufficient number of points
  n <- sum(counts.all < min.points)
  if (n == 0) {
    counts <- counts.all
    ids.fewpoints <- integer(0)

  } else {
    ids.fewpoints <- ids[c(counts.all) < min.points]
    ids <- setdiff(ids, ids.fewpoints)

    if (length(ids) == 0) stop("No flight lines with sufficient points")

    counts <- counts.all[ as.character(ids) ]
  }

  counts <- data.frame(
    flightlineID = ids,
    npoints = counts
  )


  rects <- lapply(ids, function(id) {
    r <- get_bounding_rectangle(las, classes = "all", flightlines = id)
    r$flightlineID <- id
    r
  })

  geoms <- lapply(rects, function(r) {
    v <- as.matrix(r[, c("X", "Y")])
    colnames(v) <- c("x", "y")
    sf::st_polygon(list(v))
  })

  sfdat <- sf::st_sf(flightlineID = ids, geometry = geoms)


  rects <- dplyr::bind_rows(rects) %>%
    dplyr::left_join(counts, by = "flightlineID") %>%

    dplyr::group_by(flightlineID) %>%

    dplyr::do({
      # Find side lengths and angle of longest side
      dxs <- diff(.$X[1:3])
      dys <- diff(.$Y[1:3])

      side12.len <- sqrt(dxs[1]^2 + dys[1]^2)
      side23.len <- sqrt(dxs[2]^2 + dys[2]^2)

      side.angle <- ifelse(
        side12.len > side23.len, atan2(dys[1], dxs[1]), atan2(dys[2], dxs[2])
      )

      # convert to degrees and express as angle to horizontal
      side.angle <- 180 * side.angle / pi
      if (side.angle < 0) side.angle <- side.angle + 180
      if (side.angle > 90) side.angle <- 180 - side.angle

      # Initially classify orientation
      o <- base::cut(
        side.angle,
        breaks = c(-1, angular.tol, 90 - angular.tol, 91),
        labels = c("EW", "XX", "NS"))

      # Convert to character to avoid factor level problems later
      orientation.initial <- as.character(o)

      # Side length ratio
      lens <- c(side12.len, side23.len)
      ratio.sidelen <- max(lens) / min(lens)

      data.frame(angle = side.angle, orientation.initial, ratio.sidelen,
                 stringsAsFactors = FALSE)
    }) %>%

    dplyr::ungroup()


  # GPS point time trend along X and Y dimensions for each flight line
  # and orientation inferred from this
  timetrend <- las@data %>%
    dplyr::filter(flightlineID %in% ids) %>%

    dplyr::select(X, Y, flightlineID, gpstime) %>%

    dplyr::group_by(flightlineID) %>%

    dplyr::mutate(dtime = gpstime - min(gpstime)) %>%

    dplyr::do({
      nrecs <- nrow(.)

      if (nrecs > 1000) ii <- sample.int(nrecs, 1000)
      else ii <- 1:nrecs

      dat <- .[ii, ]

      model <- lm(dtime ~ X + Y, data = dat)
      pdat <- expand.grid(X = range(dat$X), Y = range(dat$Y)) %>% dplyr::arrange(X, Y)
      p <- predict(model, newdata = pdat)
      dtX <- abs(p[1] - p[3])
      dtY <- abs(p[1] - p[2])

      otime <- ifelse(dtY > 2 * dtX, "NS", ifelse(dtX > 2 * dtY, "EW", "XX"))

      data.frame(dtX, dtY, orientation.time = otime,
                 stringsAsFactors = FALSE)
    })


  rects <- dplyr::left_join(rects, timetrend,
                            by = "flightlineID") %>%

    dplyr::mutate(orientation = dplyr::case_when(
      orientation.initial == "XX" ~ orientation.initial,
      ratio.sidelen > min.ratio ~ orientation.initial,
      TRUE ~ orientation.time
    )) %>%

    dplyr::select(flightlineID, angle, ratio.sidelen, orientation)


  sfdat <- sfdat %>%
    dplyr::left_join(rects, by = "flightlineID") %>%
    dplyr::left_join(counts, by = "flightlineID") %>%
    dplyr::mutate(ppoints = npoints / total.count)

  if (length(ids.fewpoints) > 0) {
    npoints.few <- counts.all[ as.character(ids.fewpoints) ]
    ppoints.few <- npoints.few / total.count

    suppressWarnings(
      sfdat.few <- sf::st_sf(
        flightlineID = ids.fewpoints,
        angle = NA_real_,
        ratio.sidelen = NA_real_,
        orientation = NA_character_,
        npoints = npoints.few,
        ppoints = ppoints.few,
        geometry = sf::st_sfc( sf::st_polygon() ),

        stringsAsFactors = FALSE)
    )

    # Note: using dplyr::rbind does not work here for some reason
    # (with dplyr 0.8.0.1)
    sfdat <- rbind(sfdat, sfdat.few) %>%
      dplyr::arrange(flightlineID)
  }

  sfdat
}


#' Find the minimum bounding rectangle that encloses specified points.
#'
#' When called with a single argument for the LAS tile, this function returns
#' the minimum bounding rectangle for the point cloud. The \code{classes} and
#' \code{flightlineIDs} arguments can be used to define a subset of points
#' to consider. If another kind of subset is required, the helper function
#' \code{link{.min_rectangle}} can be used directly.
#'
#' @param las A LAS object, e.g. imported using \code{\link{prepare_tile}}.
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
#' @importFrom dplyr %>%
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


#' Remove overlap between flight lines.
#'
#' This function identifies areas of overlap between pairs of flight lines,
#' places a boundary along (approximately) the middle of the overlap area, and
#' removes points from each side that belong to the flight line (mainly) on the
#' other side. Before searching for, and removing, overlaps the function checks
#' that all flight lines have consistent orientation, ie. either all are
#' north-south or all are east-west. If this is not the case a warning message
#' is issued and the original, unmodified tile is returned.
#'
#' We added this function because LIDAR data provided for New South Wales in
#' 2018 had overlap between flight lines removed for all point classes other
#' than class 5 (high vegetation).
#'
#' @param las A LAS object, e.g. imported using \code{\link{prepare_tile}}.
#'
#' @param classes Vector of one or more integer codes for point classes to consider.
#'
#' @param res Raster cell size to use when identifying overlap regions. Default
#'   value is 10 which assumes that map units are metres.
#'
#' @param buffer Width of a buffer (map units) placed around the tile to ensure
#'   that the boundaries partitioning overlap areas extend beyond all points.
#'   Default value is 100 which assumes map units are metres.
#'
#' @param ... Additional arguments to pass to \code{link{get_flightline_info}}
#'   when checking for consistent north-south or east-west orientation.
#'
#' @param min.points The minimum number of points in a flight line for it to be
#'   considered. Flight lines with fewer points are discarded. The default value
#'   (1000) is intended to exclude flight lines that only appear at the margins
#'   of the tile. Normally, such sparse flight lines will be removed when
#'   importing the tile.
#'
#' @return A modified copy of the input LAS object.
#'
#' @seealso \code{\link{get_flightline_info}}, \code{\link{check_flightlines}},
#'   \code{\link{plot_flightlines}}
#'
#' @importFrom dplyr %>%
#'
#' @export
#'
remove_flightline_overlap <- function(las,
                                      classes = 3:5, res = 10, buffer = 100,
                                      min.points = 1000,
                                      ...) {

  flines <- sort(unique(las@data$flightlineID))
  nlines <- length(flines)

  if (nlines == 1) {
    message("Tile only contains one flight line. Nothing to do.")
    return(las)
  }

  # Get flight line info (bounding rectangles etc) and check for
  # consistent orientation
  fline.dat <- get_flightline_info(las, min.points = min.points, ...)

  # Check that each flight line has sufficient points
  fline.dat <- fline.dat %>%
    dplyr::mutate(included = npoints >= min.points)

  if (!any(fline.dat$included)) {
    stop("No flight lines have sufficient points.")
  }

  # Remove points for flight lines with too few points
  ids <- fline.dat$flightlineID[ fline.dat$included ]
  ii <- las@data$flightlineID %in% ids
  las@data <- las@data[ii, ]

  fline.dat <- dplyr::filter(fline.dat, included)

  if (nrow(fline.dat) == 1) {
    message("Only one flight line with sufficient points was retained.")
    return(las)
  }


  # At least two flight lines have sufficient points...
  #
  o <- fline.dat$orientation[ fline.dat$included ]
  ok <- (o[1] %in% c("NS", "EW")) && all(o == o[1])
  if (!ok) {
    warning("Flight lines have inconsistent orientation. Cannot check overlaps.")
    return(las)
  }

  # Better variable name for code below
  orientation <- o[1]

  # Record extent of tile
  tile.xlims <- range(las@data$X) + c(-buffer, buffer)
  tile.ylims <- range(las@data$Y) + c(-buffer, buffer)

  # Indices of points to examine
  irecs <- which(las@data$Classification %in% classes)

  # Flag vector to mark points for removal
  to.remove <- logical(length(irecs))

  # Subset of point data for processing
  dat <- as.matrix( las@data[irecs, c("X", "Y", "flightlineID")] )

  # All pairwise combinations of flightlines
  fline.pairs <- utils::combn(flines, 2)


  for (i in 1:ncol(fline.pairs)) {
    fline1 <- fline.pairs[1,i]
    fline2 <- fline.pairs[2,i]

    # Is there an overlap between this pair of lines?
    indices <- match(c(fline1, fline2), fline.dat$flightlineID)
    rect1 <- fline.dat$geometry[[indices[1]]]
    rect2 <- fline.dat$geometry[[indices[2]]]
    is.overlap <- sf::st_intersects(rect1, rect2, sparse=FALSE)[1,1]

    if (is.overlap) {
      active <- dat[,3] %in% c(fline1, fline2)

      xyf <- dat[active, ]

      remove <- .do_remove_overlap(xyf, res, tile.xlims, tile.ylims, orientation)

      # Flag removal of points from the tile. The logical OR `|` is to take into
      # account that a point might have already been flagged for removal when
      # looking at a previous pair of flightlines.
      to.remove[active] <- remove | to.remove[active]
    }
  }

  # Remove flagged points from tile.
  remove.recs <- irecs[to.remove]
  keep.recs <- setdiff(1:nrow(las@data), remove.recs)
  las@data <- las@data[keep.recs, ]

  las
}


# Private helper function for remove_flightline_overlaps
#
.do_remove_overlap <- function(xyf, res, tile.xlims, tile.ylims, orientation) {
  lines <- sort(unique(xyf[,3]))
  if (length(lines) != 2) stop("Bummer: bad data passed to .do_remove_overlap")

  fline1 <- lines[1]
  fline2 <- lines[2]

  # Call C++ function to locate approximate median X for
  # Y increments of `res`

  mids <- get_overlap_midpoints(xyf, res, orientation == "NS")
  # mids <- TEST_get_overlap_midpoints(xyf, res, orientation == "NS")

  mids <- as.data.frame(mids)
  colnames(mids) <- c("x", "y", "n")

  # Check if there are many missing values for mid-point ordinates.
  # This happens when there is no or very little overlap between
  # the flight lines, or where one line only spans are small part
  # of the tile. In either case, we will leave the data as-is.
  #
  is.na.row <- apply(mids, 1, function(xs) any(is.na(xs)))
  if (sum(is.na.row) > nrow(mids) / 2) {
    # No need to remove any points
    remove <- rep(FALSE, nrow(xyf))

  } else {
    # Found an adequate number of overlap mid-points.
    # Derive a boundary by fiting a smoothing spline.
    if (orientation == "NS") {
      m <- mgcv::gam(x ~ s(y), data = mids)
      ys <- seq(tile.ylims[1], tile.ylims[2], length.out = 100)
      xs <- mgcv::predict.gam(m, newdata = data.frame(y = ys))
    } else { # "EW"
      m <- mgcv::gam(y ~ s(x), data = mids)
      xs <- seq(tile.xlims[1], tile.xlims[2], length.out = 100)
      ys <- mgcv::predict.gam(m, newdata = data.frame(x = xs))
    }

    # browser()

    # Identify points on one side of the boundary
    if (orientation == "NS") {
      pol.x <- c(xs, tile.xlims[1], tile.xlims[1], xs[1])
      pol.y <- c(ys, tile.ylims[2], tile.ylims[1], ys[1])
    } else { #EW
      pol.x <- c(xs, tile.xlims[2], tile.xlims[1], xs[1])
      pol.y <- c(ys, tile.ylims[2], tile.ylims[2], ys[1])
    }

    inside <- sp::point.in.polygon(xyf[,1], xyf[,2], pol.x, pol.y) > 0

    # Identify flight line to assign on either side of the boundary
    # The following steps try to allow for the fact that both lines
    # might be more or less equally present on one side of the
    # boundary
    tbl <- table(xyf[,3], inside)
    tbl <- apply(tbl, 2, function(vals) vals / sum(vals))
    top <- which.max(tbl)
    if (top == 2 | top == 3) {
      in.fline <- fline1
    } else {
      in.fline <- fline2
    }

    # Remove a point if it is inside the polygon (on the reference side of the
    # boundary) but does not belong to the inside flightline, OR it is outside
    # (other side of the boundary) but belongs to the inside flightline.
    remove <- inside != (xyf[,3] == in.fline)
  }

  remove
}



#' Plot point locations and flight lines
#'
#' Displays points locations for a random sample of points from a LAS object,
#' with points coloured by flight line ID. This is useful for detecting
#' overlapping flightlines.
#'
#' @param las A LAS object, e.g. imported using \code{\link{prepare_tile}}.
#'
#' @param npts Number of points to draw (default: 5000). The subset of
#'   points is selected to cover all combinations of flight line and
#'   point class.
#'
#' @param shape Point shape, specified as an integer code as for gpplot and
#'   R base plot.
#'
#' @param size Point size.
#'
#' @return A ggplot object.
#'
#' @importFrom dplyr %>%
#' @importFrom ggplot2 aes coord_equal element_blank ggplot geom_point theme
#'
#' @examples
#' \dontrun{
#' # Quick plot
#' plot_flightlines(las)
#'
#' # Facet by point classes
#' plot_flightlines(las) + facet_wrap(~ Classification)
#' }
#'
#' @export
#'
plot_flightlines <- function(las, npts = 5000, shape = 16, size = 1) {
  prop <- min(1.0, 5000 / nrow(las@data))

  if ("flightlineID" %in% colnames(las@data)) {
    dat <- las@data %>%
      as.data.frame() %>%
      dplyr::group_by(flightlineID, Classification) %>%

      dplyr::do({
        N <- nrow(.)
        n <- max(1, prop * N)
        .[sample.int(N, n), ]
      }) %>%

      dplyr::ungroup() %>%

      dplyr::mutate(flightline = factor(flightlineID))

    ggplot(data = dat, aes(x = X, y = Y)) +
      geom_point(aes(colour = flightline),
                 shape = shape, size = size) +

      coord_equal() +

      theme(axis.ticks = element_blank(),
            axis.text = element_blank())
  }
  else {
    warning("No flightlineID field found")
  }
}


#' Get scan times for a LAS tile
#'
#' The data table for an imported LAS tile (a \code{las} object) includes a
#' \code{gpstime} column which gives, for each point, scan time expressed as
#' \code{S - 1e9} where \code{S} is the number of seconds since GPS epoch time:
#' 1980-01-06 00:00:00 (GMT / UTC). This function converts the GPS time values
#' to \code{POSIXct} date-times and finds the start and end values, either for
#' the tile as a whole or for individual flightlines.
#'
#' @param las A LAS object, e.g. imported using \code{\link{prepare_tile}}.
#'
#' @param by One of 'all' (default) or 'flightline'. Case-insensitive and may
#'   be abbreviated.
#'
#' @return A data frame of start and end times.
#'
#' @importFrom dplyr %>%
#' @export
#'
get_scan_times <- function(las, by) {
  if (missing(by) || is.null(by) || by == "") {
    by <- "all"
  } else {
    by <- match.arg(tolower(by), c("all", "flightline"))
  }

  T0 <- lubridate::ymd_hms("1980-01-06 00:00:00", tz = "GMT")
  times <- T0 + lubridate::seconds(1e9 + las@data$gpstime)

  if (by == "all") {
    data.frame(time.start = min(times), time.end = max(times))

  } else {
    data.frame(flightlineID = las@data$flightlineID, times) %>%
      dplyr::group_by(flightlineID) %>%
      dplyr::summarize(time.start = min(times), time.end = max(times)) %>%
      dplyr::arrange(flightlineID)
  }
}


#' Check whether a LAS tile object has been prepared
#'
#' When a LAS tile is imported using the \code{\link{prepare_tile}} function,
#' point heights are normalized relative to ground elevation. This adds an
#' extra column (Zref) to the LAS data table.
#'
#' @param las A LAS object, e.g. imported using \code{\link{prepare_tile}}.
#'
#' @return TRUE if the tile has been prepared or FALSE otherwise.
#'
#' @export
#'
is_prepared_tile <- function(las) {
  if (!inherits(las, "LAS")) stop("Object is not a LAS tile")
  all(c("zref", "flightlineid") %in% tolower(colnames(las@data)))
}


#' Get vegetation point counts for defined strata.
#'
#' @param las A LAS object, e.g. imported using \code{\link{prepare_tile}}.
#'
#' @param strata A data frame of strata definitions with three columns:
#'   name, lower, upper.
#'
#' @param res Raster cell size.
#'
#' @param classes Integer classification codes for points to include.
#'   Default is 2 (ground) and 3:5 (vegetation).
#'
#' @return A \code{RasterStack} where each layer corresponds to a stratum
#'   (or ground) and cell values are point counts.
#'
#' @export
#'
get_stratum_counts <- function(las, strata, res = 10, classes = c(2,3,4,5)) {

  if (!is_prepared_tile(las))
    stop("Argument 'las' must be a LAS object imported with prepare_tile\n",
         "  and with point heights normalized to ground level")

  strata <- check_strata(strata)

  # Convert lower and upper strata heights to a set of cut-points
  # for use with the findInterval function.
  # This code assumes all upper heights are unique which should have
  # been checked by check_strata().
  #
  n <- nrow(strata)
  b <- unique(sort(c(strata$lower, strata$upper)))

  stratabreaks = list(
    nstrata = n,
    breaks = b,
    strata.indices = match(b, strata$upper)[-1]
  )

  xy <- lidR:::group_grid(las@data$X, las@data$Y, res)
  xy <- matrix(c(xy[[1]], xy[[2]]), ncol = 2)
  colnames(xy) <- c("X", "Y")

  r <- raster(
    xmn = min(xy[,1]) - res/2,
    xmx = max(xy[,1]) + res/2,
    ymn = min(xy[,2]) - res/2,
    ymx = max(xy[,2]) + res/2,
    res = res
  )

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
      if (nrow(dat) == 0) setValues(r, 0)
      else rasterize(dat[, 1:2], r, fun = 'count', background = 0)
    })

  rcounts <- stack(rcounts)
  names(rcounts) <- strata$name

  rcounts
}


#' Extract building points from a LAS tile
#'
#' @param las A LAS object, e.g. imported using \code{\link{prepare_tile}}.
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
#' @importFrom dplyr %>%
#'
#' @export
#'
get_building_points <- function(las) {
  if (!is_prepared_tile(las))
    stop("Argument 'las' must be a LAS object imported with prepare_tile")

  dat <- las@data %>%
    as.data.frame() %>%
    dplyr::filter(Classification == 6) %>%
    dplyr::select(X, Y, Z)

  sf::st_as_sf(dat, coords = c("X", "Y"))
}


#' Calculate strata cover estimates
#'
#' This function takes a \code{RasterStack} of point count data for the ground
#' layer and vegetation strata and calculates corresponding layer cover
#' estimates. It assumes that the first layer is for ground points and
#' subsequent layers are for vegetation strata in increasing height order.
#'
#' @param rcounts A \code{RasterStack} or \code{RasterBrick} with a layer for
#'   each vegetation stratum plus a layer for ground. Cell values are point
#'   counts.
#'
#' @return A \code{RasterStack}
#'
#' @export
#'
get_stratum_cover <- function(rcounts) {
  N <- nlayers(rcounts)
  if (N < 2) stop("Expected at least two layers: ground plus one or more strata")

  rsum <- rcounts[[1]]
  rcover <- lapply(2:N, function(i) {
    rsum <<- rsum + rcounts[[i]]
    r <- rcounts[[i]] / rsum
    r[is.nan(r)] <- 0
    r
  })

  names(rcover) <- names(rcounts)[-1]

  stack(rcover)
}


#' Metrics function for summary statistics
#'
#' @importFrom dplyr %>%
#'
#' @export
#'
layer_summary_metrics <- function(z, prob, breaks) {
  labels <- 1:(length(breaks)-1)

  zcat <- cut(
    z,
    breaks = breaks,
    labels = labels,
    right = TRUE)

  dat <- data.frame(z, zcat)

  x <- dat %>%
    dplyr::group_by(zcat) %>%
    dplyr::summarize(n = n(),
              mean = mean(z),
              median = median(z),
              lwr = highest_density_interval(z, 0.5)[1],
              upr = highest_density_interval(z, 0.5)[2])

  if (length(z) == 1) {
    bounds <- c(z, z)
  } else {
    bounds <- unname(hpdi.vec(z, prob))
  }

  list(mean = mean(z),
       median = median(z),
       lwr = bounds[1],
       upr = bounds[2])
}


#' Calculate highest density interval for a vector of values
#'
#' @export
#'
highest_density_interval <- function (x, prob = 0.95) {
  n <- length(x)
  if (n <= 1) stop("x must have more than 1 element")
  x <- sort(x)

  gap <- max(1, min(n - 1, round(n * prob)))
  init <- 1:(n - gap)

  inds <- which.min(x[init + gap] - x[init])

  out <- c(lower = x[inds], upper = x[inds + gap])
  out
}


#' Check strata lookup table and put it into a standard form
#'
#' This function is used by \code{\link{get_stratum_counts}}. It checks that:
#' \itemize{
#'   \item The object is a data frame with (at least) three columns labelled:
#'     name, lower, upper. Column order and case of names does not matter.
#'   \item All stratum upper heights are greater than lower heights.
#'   \item No two strata overlap (gaps between strata are permitted).
#' }
#' An error results if any of these conditions are not met.
#'
#' @param strata A data frame of strata definitions.
#'
#' @return A data frame in standard form with columns name, lower, upper,
#'   and rows ordered by ascending stratum heights.
#'
#' @export
#'
check_strata <- function(strata) {
  if (!inherits(strata, "data.frame"))
    stop("Strata lookup table should be a data frame")

  names <- tolower( colnames(strata) )
  if (!all( c("name", "lower", "upper") %in% names ))
    stop("Strata lookup table should have columns: name, lower, upper")

  colnames(strata) <- names
  strata <- strata[, c("name", "lower", "upper")]

  if (!all(strata$upper > strata$lower))
    stop("All upper heights should be greater than lower heights ",
         "in strata lookup table")

  o <- order(strata$lower)
  strata <- strata[o, ]

  n <- nrow(strata)
  if (!all(strata$lower[-1] >= strata$upper[-n]))
    stop("Strata lookup table has overlapping strata")

  strata
}

