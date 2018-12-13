#' Import and prepare a LAS tile for further processing
#'
#' This function imports data from a LAS file and prepares it for further
#' processing. A surface model is first fitted to ground points (class 2) by
#' Delaunay triangulation and the elevation of all points is then adjusted to be
#' relative to ground level.
#'
#' @param path Path to the LAS tile to process.
#'
#' @param drop.negative If TRUE, any points (other than ground points) whose
#'   heights are below ground level (as estimated by Delaunay interpolation of
#'   ground point heights) are discarded.
#'
#' @param fields A string containing single-letter abbreviations for the data
#'   fields to include. The default is to return all fields. See
#'   \code{\link[lidR]{readLAS}} for details of the available abbreviations.
#'
#' @param classes Point classifcations to include. Default is ground (2), vegetation
#'   (3, 4, 5) and buildings (6).
#'
#' @param flight.gap The minimum time gap (seconds) to use when assigning points
#'   to flight lines.
#'
#'
#' @return A \code{LAS} object.
#'
#' @export
#'
prepare_tile <- function(path,
                         drop.negative = TRUE,
                         fields = "*",
                         classes = c(2, 3, 4, 5, 6),
                         flight.gap = 60) {

  if (length(path) > 1) {
    warning("Presently only one path element is supported. Ignoring extra elements.")
    path <- path[1]
  }

  filtertxt <- paste("-keep_class", paste(classes, collapse = " "))

  las <- lidR::readLAS(path, select = fields, filter = filtertxt)

  # Normalize point heights relative to ground level
  lidR::lasnormalize(las, method = "delaunay")

  # Add flight line indices based on GPS times for points
  lidR::lasflightline(las, dt = flight.gap)

  # remove negative heights
  if (drop.negative) {
    ii <- las@data$Classification == 2 | las@data$Z > 0
    las@data <- las@data[ ii, ]
  }

  las
}


#' Remove overlap between flight lines.
#'
#' This function identifies areas of overlap between pairs of flight lines,
#' places a boundary along (approximately) the middle of the overlap area, and
#' removes points from each side that belong to the flight line (mainly) on the
#' other side.
#'
#' We added this function because LIDAR data provided for New South Wales in
#' 2018 has overlap between flight lines removed for all point classes other
#' than class 5 (high vegetation).
#'
#' @note The algorithm used assumes that flight lines are approximately
#'   north-south. It will fail spectacularly if this is not the case. If you
#'   encounter data where the flight lines run east-west, a (hack) work-around
#'   is to swap the X and Y ordinates before passing a LAS object to this
#'   function, then swap them back again in the returned object
#'
#' @param las A LAS object, e.g. imported using \code{\link{prepare_tile}} or
#'   \code{\link[lidR]{readLAS}}.
#'
#' @param classes Vector of one or more integer codes for point classes to consider.
#'
#' @param res Raster cell size to use when delineating overlap areas.
#'
#' @param buffer Buffer width. This is used to ensure the boundary partitioning an
#'   overlap area extends beyond all points.
#'
#' @return A modified copy of the input LAS object.
#'
#' @export
#'
remove_flightline_overlap <- function(las, classes = 5, res = 10, buffer = 100) {
  flines <- sort(unique(las@data$flightlineID))
  nlines <- length(flines)

  if (nlines > 1) {
    # Record extent of tile
    xrange <- range(las@data$X) + c(-buffer, buffer)
    yrange <- range(las@data$Y) + c(-buffer, buffer)

    # Indices of points to examine
    irecs <- which(las@data$Classification %in% classes)

    # Flag vector to mark points for removal
    to.remove <- logical(length(irecs))

    # Subset of point data for processing
    dat <- as.matrix( las@data[irecs, c("X", "Y", "flightlineID")] )


    # derive a boundary for each pair of flight lines and use it
    # to segment points in the overlap
    for (i in 1:(nlines-1)) {
      fline1 <- flines[i]
      fline2 <- flines[i+1]

      active <- dat[,3] %in% c(fline1, fline2)
      xyf <- dat[active, ]

      # Call C++ function to locate approximate median X for
      # Y increments of `res`
      mids <- get_overlap_midpoints(xyf, res)
      mids <- as.data.frame(mids)
      colnames(mids) <- c("x", "y", "n")

      # Derive a boundary by fiting a smoothing spline to the overlap mid-points
      m <- mgcv::gam(x ~ s(y), data = mids)
      ys <- seq(yrange[1], yrange[2], length.out = 100)
      xs <- mgcv::predict.gam(m, newdata = data.frame(y = ys))

      # Identify points on one side of the boundary
      pol.x <- c(xs, xrange[1], xrange[1], xs[1])
      pol.y <- c(ys, yrange[2], yrange[1], ys[1])

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

      # Flag removal of points from the tile. The logical OR `|` is to take into
      # account that a point might have already been flagged for removal when
      # looking at a previous pair of flightlines.
      to.remove[active] <- to.remove[active] | remove
    }

    # Remove flagged points from tile.
    remove.recs <- irecs[to.remove]
    keep.recs <- setdiff(1:nrow(las@data), remove.recs)
    las@data <- las@data[keep.recs, ]
  }

  las
}


#' Plot point locations and flight lines
#'
#' Displays points locations for a random sample of points from a LAS object,
#' with points coloured by flight line ID. This is useful for detecting
#' overlapping flightlines.
#'
#' @param las A LAS object, e.g. imported using \code{\link{prepare_tile}} or
#'   \code{\link[lidR]{readLAS}}.
#'
#' @param npts Number of points to sample from the LAS object (default: 5000).
#'
#' @param shape Point shape, specified as an integer code as for gpplot and
#'   R base plot.
#'
#' @param size Point size.
#'
#' @return A ggplot object.
#'
#' @importFrom dplyr %>%
#' @importFrom ggplot2 aes coord_equal ggplot geom_point
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
  if ("flightlineID" %in% colnames(las@data)) {
    dat <- las@data %>%
      as.data.frame() %>%
      dplyr::sample_n(npts) %>%
      dplyr::mutate(flightline = factor(flightlineID))

    ggplot(data = dat, aes(x = X, y = Y)) +
      geom_point(aes(colour = flightline),
                 shape = shape, size = size) +

      coord_equal()
  }
  else {
    warning("No flightlineID field found")
  }
}


#' Check whether a LAS tile object has been prepared
#'
#' When a LAS tile is imported using the \code{\link{prepare_tile}} function,
#' point heights are normalized relative to ground elevation. This adds an
#' extra column (Zref) to the LAS data table.
#'
#' @param las A LAS tile, ie. an imported object rather than a path.
#'
#' @return TRUE if the tile has been prepared or FALSE otherwise.
#'
#' @export
#'
is_prepared_tile <- function(las) {
  if (!inherits(las, "LAS")) stop("Object is not a LAS tile")
  "zref" %in% tolower(colnames(las@data))
}


#' Get vegetation point counts for defined strata.
#'
#' @param las A LAS tile imported and prepared with function
#'   \code{\link{prepare_tile}}.
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
get_strata_counts <- function(las, strata, res = 10, classes = c(2,3,4,5)) {

  if (!is_prepared_tile(las))
    stop("Argument 'las' must be a LAS object imported with prepare_tile")

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
#' @param las A LAS tile imported and prepared with function
#'   \code{\link{prepare_tile}}.
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
#' @importFrom sf st_as_sf
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

  st_as_sf(dat, coords = c("X", "Y"))
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
#' @importFrom dplyr %>% group_by
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

  x <- dat %>% group_by(zcat) %>%
    summarize(n = n(),
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
#' This function is used by \code{\link{get_strata_counts}}. It checks that:
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


