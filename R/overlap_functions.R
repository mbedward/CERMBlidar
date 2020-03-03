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
#' than class 5 (high vegetation). Note that in some cases there can be overlap
#' between bounding rectangles of flight lines returned by the function
#' \code{get_flightline_info} but no overlap for points in the class(es) being
#' considered by this function.
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
#' @param min.points The minimum number of points in the relevant point classes
#'   for a flight line to be considered. Flight lines with fewer points are
#'   ignored. If less than two flight lines have sufficient points the LAS
#'   object is returned unchanged.
#'
#' @return A modified copy of the input LAS object.
#'
#' @seealso \code{\link{get_flightline_info}}, \code{\link{check_flightlines}},
#'   \code{\link{plot_flightlines}}
#'
#' @export
#'
remove_flightline_overlap <- function(las,
                                      classes = 5, res = 10, buffer = 100,
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
  fline.dat <- get_flightline_info(las, min.points = 0, ...)

  # Check that each flight line has sufficient points in the classes being
  # considered
  x <- table(las@data[, c("flightlineID", "Classification")])[, as.character(classes), drop=FALSE]
  x <- rowSums(x) >= min.points

  fline.dat <- fline.dat %>%
    dplyr::mutate(included = x)

  if (!any(fline.dat$included)) {
    warning("No flight lines have sufficient points in relevant classes. LAS object is unchanged.")
    return(las)
  }

  # Remove any flight lines with too few points in the relevant classes
  # from further consideration
  fline.dat <- dplyr::filter(fline.dat, included)
  flines <- fline.dat$flightlineID

  if (nrow(fline.dat) == 1) {
    message("Only one flight line with sufficient points was retained.")
    return(las)
  }


  # At least two flight lines have sufficient points...
  #
  o <- fline.dat$orientation[ fline.dat$included ]
  ok <- (o[1] %in% c("NS", "EW")) && all(o == o[1])
  if (!ok) {
    stop("Flight lines have inconsistent orientation. Cannot check overlaps.")
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

      # Further check that we have two flight lines in the subset of data.
      # This is required because bounding polygons of the flight lines are based
      # on all point classes whereas here we are (usually) only considering
      # one or a few classes.
      if (length(unique(xyf[,3])) > 1) {
        remove <- .do_remove_overlap(xyf, res, tile.xlims, tile.ylims, orientation)

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



#' Remove overlap between flight lines using a raster approach
#'
#' @export
#'
remove_flightline_overlap2 <- function(las, res = 10, classes = 5) {
  r <- .nibble_flightlines(las, res = res)$nibbled
  target.id <- raster::extract(r, cbind(las@data$X, las@data$Y))

  target <- las@data$Classification %in% classes
  keep <- !target | las@data$flightlineID == target.id

  las@data <- las@data[keep, ]
  update_tile_header(las)
}


# Return the majority value from a vector. If two or more values
# are tied, choose randomly.
.majority <- function(x, ...) {
  if (all(is.na(x))) NA
  else {
    h <- sort(table(x))
    x <- as.numeric(names(h))
    ii <- which(h == max(h))
    if (length(ii) == 1) x[ii]
    else x[sample(ii, 1)]
  }
}


# Create a raster stack with a layer for each flight line and
# point counts as cell values.
.make_idcount_layers <- function(las, res = 10, classes = NULL) {
  stopifnot("flightlineID" %in% colnames(las@data))

  idmax <- max(las@data$flightlineID, na.rm = TRUE)

  f <- function(ids, ...) {
    tabulate(ids, nbins = idmax)
  }

  e <- extent(get_las_bounds(las)[1:4])

  if (is.null(classes)) classes <- unique(las@data$Classification)
  ii <- las@data$Classification %in% classes

  r <- rasterize(cbind(las@data$X[ii], las@data$Y[ii]),
                 raster(e, res = res),
                 field = las@data$flightlineID[ii],
                 fun = f,
                 background = NA)
  r
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
.nibble_flightlines <- function(las, res = 10, classes = NULL) {
  r.idcounts <- .make_idcount_layers(las, res = res, classes = classes)
  r.base <- .make_idbase_layer(r.idcounts)
  r <- r.base
  w <- matrix(1, 3, 3)

  nas <- anyNA(values(r))
  if (nas) {
    w <- matrix(1, 3, 3)
    while(nas) {
      r <- raster::focal(r, w, fun = .majority, na.rm = TRUE, pad = TRUE, NAonly = TRUE)
      nas <- anyNA(values(r))
    }
  }

  list(base = r.base, nibbled = r)
}
