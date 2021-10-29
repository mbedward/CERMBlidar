#' Correct biased point class frequency in overlap between flight lines
#'
#' This function can be used to correct for biased representation of one or more
#' point classes in areas where adjacent flight lines overlap. Such biases can
#' be caused by data post-processing where flight line overlap is removed for
#' some, but not all, point classes. Why anyone would want to do that is a
#' complete mystery, but it was commonly seen in airborne LiDAR data for New
#' South Wales, Australia, prior to 2020 where overlap was removed for all but
#' class 5 (high vegetation) points resulting in incorrect vegetation cover
#' estimates. If the overlap is either left for all point classes, or removed
#' for all classes, then vegetation cover and other ratio-based metrics will be
#' unbiased.
#'
#' This function identifies areas of overlap between adjacent flight lines based
#' on the specified point classes; determines a mid-line boundary between the
#' flight lines within each area; and removes points from the data such that the
#' remaining points from each flight line are strictly on one side of the
#' boundary. Two algorithms are available: the first ('spline') identifies a
#' vector boundary in continuous space; the second ('nibble') uses a raster
#' approach that is more robust in the face of inconsistent flight line shapes
#' and orientations, but can result in minor image artefacts in the form of
#' reduced point density. By default, the spline boundary algorithm is tried
#' first and, only if this fails, the nibble algorithm is applied. See details
#' below. \strong{Caution:} You should only use this function after checking
#' that there actually is a bias in point class representation within any
#' overlap areas using the function \code{check_overlap_bias}, and be
#' careful to specify the correct classes to consider via the \code{classes}
#' argument. If no bias exists in the data, running this function for a subset
#' of point classes will create one!
#'
#' The spline boundary algorithm begins by identifying areas of overlap in
#' raster space, then fitting a spline vector boundary along the middle of each
#' area, and removes points from each side that belong to the flight line
#' (mainly) on the other side. Before searching for, and removing, overlaps the
#' function checks that all flight lines have consistent orientation: either all
#' north-south or all are east-west, allowing for some leeway. If this is not
#' the case, the algorithm will not proceed.
#'
#' The nibble algorithm also begins by identifying areas of overlap in raster
#' space. Next it creates a raster where cells outside overlap areas are set to
#' the integer flight line ID, while cells in overlap areas are set to missing
#' values. A 'nibble' process (inspired by the ArcGIS raster algorithm of the
#' same name) is then iteratively applied to replace the missing values with the
#' majority value of neighbouring cells. This tends to produce a similar
#' partitioning to the vector boundary algorithm for overlap areas that span the
#' image, but also deals with irregular flight line orientations and shapes. The
#' disadvantage of this algorithm is that it can leave narrow trails of reduced
#' point density across the image if any pairs of flight lines abut rather than
#' overlap. However, this effect can be minimized by working at a finer raster
#' resolution.
#'
#' @param las A LAS object, e.g. imported using \code{prepare_tile}.
#'
#' @param algorithms A vector of algorithm names. Valid options are 'spline' and
#'   'nibble'. This argument can be used to restrict processing to just one
#'   algorithm. If both are specified (the default) the spline algorithm will
#'   always be tried first as it tends to produce better results but can fail if
#'   flight line orientation is inconsistent.
#'
#' @param classes Vector of one or more integer codes for point classes to
#'   consider. Only points belonging to these classes will be removed.
#'
#' @param res.spline Raster cell size to use when identifying flight line
#'   overlap areas using the spline algorithm. The default (10) assumes map
#'   units are metres.
#'
#' @param res.nibble Raster cell size to use when identifying and removing
#'   overlap areas using the nibble algorithm. The default is half of the
#'   \code{res.spline} value if specified, otherwise 5 (which assumes map units
#'   are metres).
#'
#' @param buffer Width of a buffer (map units) placed around the tile to ensure
#'   that, for the 'spline' algorithm, the boundaries partitioning overlap areas
#'   extend beyond all points. Default value is 100 which assumes map units are
#'   metres.
#'
#' @param ... Additional arguments to pass to \code{link{get_flightline_info}}
#'   when checking for consistent north-south or east-west orientation.
#'
#' @param min.points The minimum number of points in the relevant point classes
#'   for a flight line to be considered. Flight lines with fewer points are
#'   ignored. If less than two flight lines have sufficient points the LAS
#'   object is returned unchanged.
#'
#' @return A modified copy of the input LAS object.
#'
#' @seealso \code{\link{get_flightline_info}}, \code{\link{check_flightlines}},
#'   \code{\link{get_flightline_polygons}}, \code{\link{plot_flightline_points}}
#'
#' @export
#'
remove_flightline_bias <- function(las,
                                   algorithms = c("spline", "nibble"),
                                   classes = 5,
                                   res.spline = 10,
                                   res.nibble = res.spline / 2,
                                   buffer = 100,
                                   min.points = 1000,
                                   ...) {

  flines <- sort(unique(las@data$flightlineID))
  nlines <- length(flines)

  if (nlines == 1) {
    message("Tile only contains one flight line. Nothing to do.")
    return(las)
  }

  algorithms <- match.arg(tolower(algorithms),
                          choices = c("spline", "nibble"),
                          several.ok = TRUE)

  done <- FALSE
  if ("spline" %in% algorithms) {
    message("Trying spline boundary method...")
    res <- .do_spline_remove_bias(las,
                                  classes = classes,
                                  res = res.spline,
                                  buffer = buffer,
                                  min.points = min.points,
                                  ...)

    if (res$success) {
      las <- res$las
      done <- TRUE
    }
  }

  if (!done && "nibble" %in% algorithms) {
    message("Trying raster nibble method...")
    res <- .do_nibble_remove_bias(las,
                                  classes = classes,
                                  res = res.nibble,
                                  min.points = min.points)

    if (res$success) {
      las <- res$las
      done <- TRUE
    }
  }

  if (done) las
  else stop("Flight line overlap removal failed")
}


# Apply the spline overlap removal algorithm. Return a list with
# elements 'success' (logical) and 'las'.
#
.do_spline_remove_bias <- function(las,
                                   classes,
                                   res,
                                   buffer,
                                   min.points,
                                   ...) {

  # Get flight line info (bounding rectangles etc) and check for
  # consistent orientation
  fline.dat <- get_flightline_info(las, min.points = 0, ...)

  # Check that each flight line has sufficient points in the classes being
  # considered
  x <- table(las@data[, c("flightlineID", "Classification")])[, as.character(classes), drop=FALSE]
  x <- rowSums(x) >= min.points

  fline.dat <- fline.dat %>%
    dplyr::mutate(included = x)

  if (!any(fline.dat$included)) {
    message("No flight lines have sufficient points in relevant classes. LAS object is unchanged.")
    return(list(success = TRUE, las = las))
  }

  # Remove any flight lines with too few points in the relevant classes
  # from further consideration
  fline.dat <- dplyr::filter(fline.dat, included)
  flines <- fline.dat$flightlineID

  if (nrow(fline.dat) == 1) {
    message("Only one flight line with sufficient points was retained.")
    return(list(success = TRUE, las = las))
  }


  # At least two flight lines have sufficient points...
  #
  o <- fline.dat$orientation[ fline.dat$included ]
  ok <- (o[1] %in% c("NS", "EW")) && all(o == o[1])
  if (!ok) {
    message("Flight lines have inconsistent orientation.")
    return(list(success = FALSE, las = las))
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

      # Further check that we have two flight lines in the subset of data.
      # This is required because bounding polygons of the flight lines are based
      # on all point classes whereas here we are (usually) only considering
      # one or a few classes.
      if (length(unique(xyf[,3])) > 1) {
        remove <- .do_spline_remove_bias_worker(xyf,
                                                res,
                                                tile.xlims,
                                                tile.ylims,
                                                orientation)

        # Flag removal of points from the tile. The logical OR `|` is to
        # take into account that a point might have already been flagged
        # for removal when looking at a previous pair of flightlines.
        to.remove[active] <- remove | to.remove[active]
      }
    }
  }

  # Remove flagged points from tile.
  if (any(to.remove)) {
    remove.recs <- irecs[to.remove]
    keep.recs <- setdiff(1:nrow(las@data), remove.recs)
    las@data <- las@data[keep.recs, ]
    las <- update_tile_header(las)
  }

  list(success = TRUE, las = las)
}


# Private helper function for do_spline_remove_bias
#
.do_spline_remove_bias_worker <- function(xyf,
                                          res,
                                          tile.xlims,
                                          tile.ylims,
                                          orientation) {
  lines <- sort(unique(xyf[,3]))
  if (length(lines) != 2) stop("Bummer: bad data passed to .do_remove_bias")

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



# Apply the nibble overlap removal algorithm. Return a list with
# elements 'success' (logical) and 'las'.
#
.do_nibble_remove_bias <- function(las,
                                   classes,
                                   res,
                                   min.points) {

  # Only consider flight lines with sufficient points in the classes being
  # considered. Any with less points are ignored (ie. their points are retained).
  ii <- las@data$Classification %in% classes
  x <- table(las@data[ii, "flightlineID"])
  included.ids <- as.integer( names(x)[x > min.points] )

  if (length(included.ids) == 0) {
    message("No flight lines have sufficient points in relevant classes. LAS object is unchanged.")
    return(list(success = TRUE, las = las))
  }

  keep <- rep(TRUE, nrow(las@data))

  x <- .nibble_flightlines(las, res = res, flightline.ids = included.ids)
  target.id <- raster::extract(x$nibbled, cbind(las@data$X, las@data$Y))

  target <- las@data$Classification %in% classes &
    las@data$flightlineID %in% included.ids

  keep <- !target | las@data$flightlineID == target.id

  las@data <- las@data[keep, ]
  las <- update_tile_header(las)

  list(success = TRUE, las = las)
}


# Return the majority value from a vector of integer ID values.
# If two or more values are tied, choose randomly. Values of
# zero are never selected.
.majority_id <- function(x, ...) {
  if (all(is.na(x))) NA
  else {
    x <- na.omit(x)
    if (all(x == 0)) 0  # no neighbours have any points
    else {
      x <- x[x > 0]
      h <- sort(table(x))
      vals <- as.numeric(names(h))
      ii <- which(h == max(h))
      if (length(ii) == 1) vals[ii]
      else vals[sample(ii, 1)]
    }
  }
}


# Derive a starting raster of non-overlapping flightline IDs
# from a raster stack of point counts. Cells in overlap areas
# are given NA values.
.make_idbase_layer <- function(r.ids) {
  calc(r.ids, fun = function(x, ...) {
    b <- x > 0
    n <- sum(b, na.rm = TRUE)
    if (n == 1) which(b)
    else if (n > 1) NA
    else 0
  })
}


# Create a raster layer of non-overlapping flight lines. Areas of overlap
# are partitioned using a raster nibble approach.
.nibble_flightlines <- function(las, res = 10, flightline.ids = NULL) {
  stopifnot("flightlineID" %in% colnames(las@data))

  if (is.null(flightline.ids)) flightline.ids <- sort(unique(las@data$flightlineID))

  r.idcounts <- .make_idcount_layers(las,
                                     res = res,
                                     flightline.ids = flightline.ids)

  r.base <- .make_idbase_layer(r.idcounts)
  r <- r.base
  w <- matrix(1, 3, 3)

  nas <- anyNA(values(r))
  if (nas) {
    w <- matrix(1, 3, 3)
    while(nas) {
      r <- raster::focal(r, w, fun = .majority_id, na.rm = TRUE, pad = TRUE, NAonly = TRUE)
      nas <- anyNA(values(r))
    }
  }

  list(base = r.base, nibbled = r)
}


# Create a raster stack with a layer for each flight line and
# point counts as cell values.
.make_idcount_layers <- function(las, res, flightline.ids = NULL) {
  stopifnot("flightlineID" %in% colnames(las@data))

  if (is.null(flightline.ids)) {
    flightline.ids <- sort(unique(las@data$flightlineID))
  } else {
    flightline.ids <- sort(unique(flightline.ids))
  }

  ii <- las@data$flightlineID %in% flightline.ids

  if (!any(ii)) stop("No points in specified flight lines")

  dat <- las@data[ii, c("X", "Y", "flightlineID")]
  id.max <- max(flightline.ids, na.rm = TRUE)
  id.ind <- which((1:id.max) %in% flightline.ids)

  f <- function(ids, ...) {
    tabulate(ids, nbins = id.max)[id.ind]
  }

  e <- extent(get_las_bounds(las)[1:4])


  r <- rasterize(dat[, c("X", "Y")],
                 raster(e, res = res),
                 field = dat$flightlineID,
                 fun = f,
                 background = 0)

  names(r) <- paste("id", flightline.ids, sep = ".")

  r
}


# Derive a starting raster of non-overlapping flightline IDs
# from a raster stack of point counts. ID values are extracted
# from layer names. Cells in overlap areas are assigned NA values.
# Cells with zero point counts in all ID layers are assigned zero.
#
.make_idbase_layer <- function(r.ids) {
  flightline.ids <- as.integer(stringr::str_extract(names(r.ids), "\\d+"))

  calc(r.ids, fun = function(x, ...) {
    b <- x > 0
    n <- sum(b, na.rm = TRUE)
    if (n == 1) flightline.ids[b]
    else if (n > 1) NA
    else 0
  })
}
