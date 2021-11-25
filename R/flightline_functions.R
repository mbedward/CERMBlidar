#' Plot a sample of points using colours to indicate flight lines
#'
#' Displays a random sample of points from a LAS object, coloured by flight line
#' ID and facetted by point class. This is useful, for example, when diagnosing
#' problems with overlapping flightlines or uneven point density.
#'
#' @param las A LAS object, e.g. imported using \code{prepare_tile}.
#'
#' @param classes An integer vector specifying the point classes to display. May
#'   include 2 (ground), 3 (low veg), 4 (mid veg), 5 (high veg), 35 (all veg -
#'   the three vegetation classes combined), 6 (buildings) and 9 (water). The
#'   default is to display ground (2) and all vegetation combined (35).
#'
#' @param npts Number of points to sample from the LAS data for drawing.
#'   Sampling is done proportionally across combinations of point classes and
#'   flight lines, and the resulting number of rows in the \code{data} element
#'   of the returned ggplot object can vary slightly. If displaying many point
#'   classes, increasing this value from the default (5000) and decreasing the
#'   point size will produce nicer plots.
#'
#' @param shape Point shape, specified as an integer code as for gpplot and
#'   R base plot.
#'
#' @param size Point size.
#'
#' @return A ggplot object.
#'
#' @importFrom ggplot2 aes as_labeller coord_equal element_blank facet_wrap ggplot geom_point theme
#'
#' @examples
#' \dontrun{
#' # Default plot of points for ground (class 2) and combined
#' # vegetation (classes 3-5).
#' plot_flightline_points(las)
#'
#' # Display the vegetation classes both separately and combined
#' plot_flightline_points(las, classes = c(3:5, 35))
#' }
#'
#' @seealso \code{\link{get_flightline_polygons}}
#'
#' @export
#'
plot_flightline_points <- function(las,
                                   classes = c(2, 35),
                                   npts = 5000,
                                   shape = 16,
                                   size = 1) {

  if (!("flightlineID" %in% colnames(las@data))) {
    stop("No flightlineID field found")
  }

  prop <- min(1.0, npts / nrow(las@data))

  ClassLookup <- data.frame(
    code = c(2, 3, 4, 5, 35, 6, 9),
    label = c('ground (2)',
              'low veg (3)', 'mid veg (4)', 'high veg (5)', 'all veg (3-5)',
              'building (6)', 'water (9)'),
    stringsAsFactors = FALSE
  )

  if (length(classes) == 0) stop("One or more point classes should be specified")

  classes <- unique(classes)
  ok <- classes %in% ClassLookup$code
  if (!all(ok)) {
    stop("Unrecognized class(es): ", classes[!ok], "\n",
         "Valid options: ", paste(ClassLookup$code, collapse = ", "))
  }

  # Get data for required classes, allowing for the display of veg
  # classes separately and/or combined
  dat <- lapply(classes, function(cl) {
    if (cl == 35) {
      cldat <- las@data %>%
        as.data.frame() %>%
        dplyr::filter(Classification %in% 3:5) %>%
        dplyr::mutate(Classification = 35)
    } else {
      cldat <- las@data %>%
        as.data.frame() %>%
        dplyr::filter(Classification == cl)
    }

    cldat
  })

  dat <- dplyr::bind_rows(dat)

  dat <- dat %>%
    dplyr::group_by(flightlineID, Classification) %>%

    dplyr::do({
      N <- nrow(.)
      ns <- max(1, prop * N)
      ii <- sample.int(N, ns)
      .[ii, ]
    }) %>%

    dplyr::ungroup() %>%

    dplyr::mutate(flightline = factor(flightlineID))


  gg <- ggplot(data = dat, aes(x = X, y = Y)) +
    geom_point(aes(colour = flightline),
               shape = shape, size = size) +

    coord_equal() +

    theme(axis.ticks = element_blank(),
          axis.text = element_blank())

  ii <- ClassLookup$code %in% classes
  labels <- ClassLookup$label[ii]
  names(labels) <- ClassLookup$code[ii]
  gg <- gg + facet_wrap(~ Classification, labeller = as_labeller(labels))

  gg
}


#' Create convex polygons for flight lines
#'
#' Given an input \code{LAS} object, this function returns an \code{sf} spatial
#' data frame with a convex polygon for each flight line or, optionally, each
#' combination of flight line and point class. The default is to only consider
#' the ground class (2) and the vegetation classes (3-5). Drawing the resulting
#' polygons grouped by point class (see example below) can help in checking for
#' biased representation of point classes in overlap areas.
#'
#' @param las A LAS object, e.g. imported using \code{prepare_tile}.
#'
#' @param classes Vector of one or more integer codes for point classes to
#'   consider. Default is \code{2:5} (ground and vegetation).
#'
#' @param group_classes If \code{TRUE}, separate polygons will be created for
#'   each combination of flight line and point class. If \code{FALSE} (default),
#'   a single polygon is created for each flight line.
#'
#' @return An \code{sf} data frame with columns \code{flightlineID},
#'   \code{Classification} (if \code{group_classes} is \code{TRUE}), and
#'   \code{group_classes}. The coordinate reference system is set to that of the
#'   input \code{LAS} object, or to its horizontal component if the input object
#'   has a compound (horizontal plus vertical) CRS.
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#'
#' fl_polys <- get_flightline_polygons(las, group_classes = TRUE)
#'
#' # Draw flight line polygons separately for each point class
#' ggplot(data = fl_polys) +
#'   geom_sf(aes(fill = factor(flightlineID)), alpha = 0.4) +
#'   facet_wrap(~ Classification)
#' }
#'
#' @export
#'
get_flightline_polygons <- function(las,
                                    classes = 2:5,
                                    group_classes = FALSE) {

  stopifnot("flightlineID" %in% colnames(las@data))

  if (is.null(classes)) {
    classes <- 2:5
  }

  lineIDs <- sort(unique(las@data$flightlineID))

  if (group_classes) {
    lookup <- expand.grid(id = lineIDs, cl = sort(classes))

    dats <- lapply(1:nrow(lookup), function(irow) {
      #browser()
      id <- lookup$id[irow]
      cl <- lookup$cl[irow]

      xy <- las@data[, c("X", "Y", "flightlineID", "Classification")] %>%
        as.data.frame() %>%
        dplyr::filter(flightlineID == id & Classification == cl) %>%
        dplyr::select(X, Y) %>%
        as.matrix()

      if (nrow(xy) > 0) {
        h <- grDevices::chull(xy)
        v <- xy[c(h, h[1]), ]

        poly <- sf::st_polygon(list(v))
        geom <- sf::st_sfc(poly, crs = get_horizontal_crs(las))

        sf::st_sf(flightlineID = id,
                  Classification = cl,
                  geometry = geom)
      }
    })

  } else { # not grouping classes
    dats <- lapply(lineIDs, function(id) {
      #browser()
      xy <- las@data[, c("X", "Y", "flightlineID")] %>%
        as.data.frame() %>%
        dplyr::filter(flightlineID == id) %>%
        dplyr::select(X, Y) %>%
        as.matrix()

      if (nrow(xy) > 0) {
        h <- grDevices::chull(xy)
        v <- xy[c(h, h[1]), ]

        poly <- sf::st_polygon(list(v))
        geom <- sf::st_sfc(poly, crs = get_horizontal_crs(las))

        sf::st_sf(flightlineID = id, geometry = geom)
      }

    })
  }

  do.call(rbind, dats)
}


#' Identifies overlapping and non-overlapping parts of flight line polygons.
#'
#' Given a set of flight line polygons created with function
#' \code{get_flightline_polygons}, this function creates a new set of
#' polygons representing the overlapping and non-overlapping parts. Note: at
#' present, point classes are not considered separately. If the input \code{sf}
#' data frame of polygons includes a \code{Classification} column, the polygons
#' for each class will be merged prior to identifying overlap and non-overlap
#' parts.
#'
#' @param polys An \code{sf} data frame of polygons as returned by function
#'   \code{get_flightline_polygons}.
#'
#' @return An \code{sf} data frame of polygons for the separate and overlapping
#'   parts of flight lines, with columns: \code{overlap} (character label; either
#'   'overlapping' or 'separate'), \code{flightlineIDs} (character label; single integer
#'   for separate parts, or colon-separated integers for overlapping parts, e.g. \code{'2:4'});
#'   \code{area} polygon area in map units; \code{geometry}.
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#'
#' fl_polys <- get_flightline_polys(las)
#' ov_polys <- get_flightline_overlaps(fl_polys)
#'
#' # Calculate polygon centroids to position labels
#' centroids <- sf::st_centroid(ov_polys)
#'
#' # Draw polygons for overlapping and non-overlapping areas
#' ggplot(data = ov_polys) +
#'   geom_sf() +
#'   geom_sf_text(data = centroids, aes(label = flightlineIDs)) +
#'   coord_sf(datum = st_crs(ov_polys)) +
#'   labs(x = "", y = "") +
#'   facet_wrap(~overlap)
#' }
#'
#' @export
#'
get_flightline_overlaps <- function(polys) {
  stopifnot("flightlineID" %in% colnames(polys))

  # If there are separate polygons for point classes, merge them
  if ("Classification" %in% colnames(polys)) {
    polys <- polys %>%
      dplyr::group_by(flightlineID) %>%
      dplyr::summarize()
  }

  # Helper function to create polygon labels for separate and overlapping
  # parts of flightlines
  fn <- function(polys, origins) {
    i <- unlist(origins)
    paste(polys$flightlineID[i], collapse = ":")
  }

  # Self-intersect the flight line polygons, then label the resulting polygons
  # as separate or overlapping and add their respective flight line IDs
  ov_polys <- sf::st_intersection(polys) %>%
    dplyr::arrange(flightlineID) %>%
    dplyr::mutate(area = as.numeric(st_area(.)),
                  overlap = factor(n.overlaps > 1, levels = c(FALSE, TRUE), labels = c("separate", "overlapping"))) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(flightlineIDs = fn(polys, origins)) %>%
    dplyr::ungroup() %>%

    dplyr::select(overlap, flightlineIDs, area, geometry)

  ov_polys
}


#' Performs a simple check for point class bias within flightline overlaps
#'
#' This function examines ground (class 2) and vegetation (classes 3-5) points
#' and checks whether the ratio of vegetation to ground points differs
#' substantially between overlapping and non-overlapping areas of flight lines.
#'
#' The check performed by this function is a simple heuristic rather than being
#' statistically rigorous, but it seems to work. For each of overlapping and
#' non-overlapping flight line areas, the number of points for each class is
#' determined. Next, the ratio of vegetation points to ground points is
#' calculated for each of the vegetation classes. Finally, these ratios are
#' compared between overlap and non-overlap areas. If any vegetation class ratio
#' in the overlap area is more than \code{bias_threshold} times greater than the
#' corresponding value in the non-overlap area, the check returns \code{TRUE}.
#'
#' @param las A LAS object, e.g. imported using \code{prepare_tile}.
#'
#' @param ov_polys An \code{sf} data frame of polygons as returned by the
#'   function \code{get_flightline_overlaps}.
#'
#' @param n_sample_points The number of points to sample from the LiDAR point
#'   cloud.
#'
#' @param bias_threshold The threshold value, above which bias is considered to
#'   be present. Default value is 2. See Details for more explanation.
#'
#' @return A logical value, where \code{TRUE} indicates that bias has been
#'   detected. An attributes list is attached to the value with elements
#'   \code{nsample_points} and \code{ratio_data}. The latter element is a
#'   data frame with a record for each point class, and columns for the
#'   number of points and ratio of vegetation class to ground class points
#'   within each of overlapping and non-overlapping areas.
#'
#' @examples
#' \dontrun{
#' fl_polys <- get_flightline_polys(las)
#' ov_polys <- get_flightline_overlaps(fl_polys)
#' res <- check_overlap_bias(las, ov_polys)
#'
#' if (res) {
#'   cat("Flight line overlap bias is present \n")
#'
#'   # Get the table of results to determine which point class(es)
#'   # have biased representation in overlap areas
#'   dat <- attr(res, "ratio_data")
#'   print(dat)
#' }
#' }
#'
#' @seealso \code{\link{remove_flightline_bias}}
#'
#' @export
#'
check_overlap_bias <- function(las,
                               ov_polys,
                               nsample_points = 1e5,
                               bias_threshold = 2.0) {

  # Create an sf point data layer for sample points
  ii <- which(las$Classification %in% 2:5)
  ii <- sample(ii, size = min(nsample_points, length(ii)))

  pts <- las@data[ii, c("X", "Y", "Classification")] %>%
    sf::st_as_sf(coords = c("X", "Y"), crs = get_horizontal_crs(las))


  # Helper function to calculate the ratio of veg class points
  # to ground points
  fn_ratio <- function(Classification, npoints) {
    iground <- which(Classification == 2)
    ratio <- npoints / npoints[iground]

    ratio
  }


  x <- sf::st_join(ov_polys, pts, join = sf::st_contains) %>%
    sf::st_drop_geometry() %>%

    # Any polygons that had no points within them will be present
    # as a record with NA for Classification. These should only be tiny
    # areas so we just drop them.
    dplyr::filter(!is.na(Classification)) %>%

    dplyr::group_by(overlap, Classification) %>%
    dplyr::summarize(npoints = n()) %>%

    dplyr::group_by(overlap) %>%
    dplyr::mutate(ratio = fn_ratio(Classification, npoints)) %>%
    dplyr::ungroup() %>%

    tidyr::pivot_wider(names_from = overlap, values_from = c(npoints, ratio))

  # Check if the ratios of veg class to ground points differ
  # substantially between overlap and separate areas
  is_bias <- any(x$ratio_overlapping / x$ratio_separate > bias_threshold)

  attributes(is_bias) <- list(nsample_points = nsample_points,
                              ratio_data = x)

  is_bias
}


#' Get summary information for flight line extents
#'
#' This function determines the minimum bounding rectangle for points in each
#' flight line and, based on this, the orientation of the flight line: one of
#' 'NS' (north - south), 'EW' (east - west), or XX (indeterminate). Note: this
#' function pre-dates the similar function \code{get_flightline_polygons} and
#' the two might be merged in a future version of the package.
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
#' @param las A LAS object, e.g. imported using \code{prepare_tile}.
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
#' @seealso \code{\link{check_flightline_orientation}}
#'
#' @export
#'
get_flightline_info <- function(las,
                                angular.tol = 25,
                                min.ratio = 1.2,
                                min.points = 1000) {

  stopifnot(angular.tol > 0.0, angular.tol <= 30)
  stopifnot(min.ratio > 1.1)

  if (is_empty_tile(las)) {
    warning("Point cloud is empty in tile")

    out <- sf::st_sf(
      flightlineID = NA_integer_,
      angle = NA_real_,
      ratio.sidelen = NA_real_,
      orientation = "XX",
      npoints = 0,
      ppoints = 0,
      geometry = sf::st_sfc( sf::st_polygon() )
    )

    return(out)
  }

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

    if (length(ids) == 0) {
      warning("No flight lines with sufficient points")
      return(NULL)
    }

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


  # Helper function to determine flight line direction in
  # the dplyr pipeline below. Also makes it easier to run
  # in a debugging session.
  fn_get_timetrend <- function(dat) {
    nrecs <- nrow(dat)

    if (nrecs > 1000) ii <- sample.int(nrecs, 1000)
    else ii <- 1:nrecs

    dat <- dat[ii, ]

    model <- lm(dtime ~ X + Y, data = dat)
    pdat <- expand.grid(X = range(dat$X), Y = range(dat$Y)) %>% dplyr::arrange(X, Y)
    p <- predict(model, newdata = pdat)
    dtX <- abs(p[1] - p[3])
    dtY <- abs(p[1] - p[2])

    otime <- ifelse(dtY > 2 * dtX, "NS", ifelse(dtX > 2 * dtY, "EW", "XX"))

    data.frame(dtX, dtY, orientation.time = otime,
               stringsAsFactors = FALSE)
  }


  # GPS point time trend along X and Y dimensions for each flight line
  # and orientation inferred from this
  timetrend <- las@data %>%
    dplyr::filter(flightlineID %in% ids) %>%

    dplyr::select(X, Y, flightlineID, gpstime) %>%

    dplyr::group_by(flightlineID) %>%

    dplyr::mutate(dtime = gpstime - min(gpstime)) %>%

    dplyr::do(fn_get_timetrend(.))


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


#' Check that the orientations of all flight lines match
#'
#' This is a convenience function that calls \code{get_flightline_info}
#' with default arguments and checks that all flight lines have the same
#' orientation, either 'EW' or 'NS'.
#'
#' @param las A LAS object, e.g. imported using \code{prepare_tile}.
#'
#' @param ... Additional arguments to pass to \code{get_flightline_info}
#'   including \code{angular.tol}.
#'
#' @return \code{TRUE} if flight lines are consistent; \code{FALSE} otherwise.
#'
#' @seealso \code{\link{get_flightline_info}}
#'
#' @export
#'
check_flightline_orientation <- function(las, ...) {
  dat <- get_flightline_info(las, ...)
  o <- na.omit(dat$orientation)

  if (length(o) == 0) FALSE
  else (o[1] %in% c("EW", "NS")) && all(o == o[1])
}

