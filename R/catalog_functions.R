#' Create a catalog of LAS tiles from HTML metadata
#'
#' This function recurses through a directory tree containing (usually) LAS
#' tiles and associated HTML metadata documents, and compiles a summary in the
#' form of a spatial data frame (class \code{sf} object). The geometry column of
#' the spatial data frame consists of rectangular polygons for tile extents.
#' Presently, the function assumes that HTML metadata is in the format used by
#' Spatial Services, New South Wales.
#'
#' @param path Top of the directory tree containing tiles to catalog.
#'
#' @param dirs One or more path elements to subset the directory tree. The
#'   default is "LAS" to ignore, for example, "DEM" directories present in data
#'   provided by Spatial Services, New South Wales. Set to NULL or an empty
#'   string if not required.
#'
#' @param exts One or more file extensions to identify metadata files.
#'
#' @param extent.crs A coordinate reference system specifier (e.g. integer EPSG
#'   code) to use for LAS tile extent polygons. The default is EPSG:4326
#'   (WGS84). Set to \code{NULL} or \code{sf::NA_crs_} if you do not want to
#'   apply a single CRS to all extent polygons.
#'
#' @return A spatial data frame containing summary metadata and tile extent
#'   polygons.
#'
#' @seealso \code{\link{read_html_metadata}}
#'
#' @importFrom sf st_polygon st_sfc st_sf
#'
#' @examples
#' \dontrun{
#' # Compile metadata and save as an ESRI shapefile
#' meta <- compile_metadata("D:/mydata")
#' sf::st_write(meta, "metadata.shp", delete_layer = TRUE)
#' }
#'
#' @export
#'
compile_metadata <- function(path, dirs = "LAS", exts = "html", extent.crs = 4326) {

  path <- normalizePath(path, winslash = "/")

  ## Build list of paths for tiles

  if (any(.str_length(exts) > 0)) {
    file.pattern <- sprintf("(%s)$", paste(exts, collapse="|"))
  } else {
    file.pattern <- NULL
  }

  files <- dir(path, recursive = TRUE, full.names = TRUE, pattern = file.pattern)

  if (length(dirs) > 0 && .str_length(dirs) > 0) {
    files <- stringr::str_subset(files, paste0(dirs, "(\\\\|/)"))
  }


  # Read and compile metadata from HTML documents
  metadata <- lapply(files, read_html_metadata)

  # Vertex indices to use with the vector of min/max longitude and latitude
  # values when constructing tile extent polygons. We follow the convention of
  # counter-clockwise vertices for exterior ring. These indices assume that
  # the vector is ordered as: lonmin, lonmax, latmin, latmax.
  #
  vi <- c(1,3, 2,3, 2,4, 1,4, 1,3)

  res <- lapply(metadata, function(m) {
    lonlat <- c(m$bounds.lon, m$bounds.lat)

    dat <- data.frame(
      filename = m$filename,
      locality = m$locality,
      ahdnum = m$ahdnum,
      startdate = as.Date(m$startdate),
      enddate = as.Date(m$enddate),
      lidnum = m$lidnum,
      classnum = m$classnum,

      lonmin = min(m$bounds.lon),
      lonmax = max(m$bounds.lon),
      latmin = min(m$bounds.lat),
      latmax = max(m$bounds.lat),

      epsg = m$epsg,

      stringsAsFactors = FALSE)

    g <- st_polygon(list(matrix(lonlat[vi], ncol=2, byrow = TRUE)))

    list(data = dat, poly = g)
  })

  # Create combined data frame
  dat <- lapply(res, function(r) r$data) %>% dplyr::bind_rows()

  # Create corresonding sf geomtry list of extent polygons

  geometry <- lapply(res, function(r) r$poly) %>% st_sfc(crs = extent.crs)

  # Return as spatial data frame
  st_sf(dat, geometry)
}


#' Read metadata for a LAS tile
#'
#' This function reads metadata from an HTML document in the format used by
#' Spatial Services, New South Wales (Australia). It also extracts attributes
#' from the document file name. It is used by function \code{compile_metadata}
#' but can also be called directly.
#'
#' @param path.html The full path to an HTML document containing tile metadata.
#'
#' @return A named list with the following elements:
#'   \describe{
#'     \item{filename}{Metadata file name.}
#'     \item{locality}{Place name (from file name).}
#'     \item{lidnum}{LID or PHO number, e.g. 'PHO3' (from file name).}
#'     \item{classnum}{Class number, e.g. 'C0' (from filename).}
#'     \item{ahdnum}{AHD number, e.g. 6386236 (from filename).}
#'     \item{startdate}{Start date of data capture (from doc).}
#'     \item{enddate}{End date of data capture (from doc).}
#'     \item{bounds.lon}{Vector of min and max longitude (from doc).}
#'     \item{bound.lat}{Vector of min and max latitude (from doc).}
#'     \item{epsg}{EPSG code for tile (from doc). To be used if converting to UTM coordinates.}
#'     \item{datum}{Horizontal datum, e.g. GDA94 (from doc).}
#'   }
#'
#' @seealso \code{\link{compile_metadata}}
#'
#' @importFrom stringr str_extract str_replace_all str_subset
#'
#' @export
#'
read_html_metadata <- function(path.html) {
  if (length(path.html) > 1)
    stop("The path.html argument is presently limited to a single file")

  if (!requireNamespace("rvest", quietly = TRUE))
    stop("You need to install the rvest package to use this function")

  ri <- function(pattern) stringr::regex(pattern, ignore_case = TRUE)

  # parse file name
  fname <- str_extract(path.html, "[^[\\\\\\/]]+$")

  locality <- fname %>% str_extract("^[A-Za-z\\s]+")

  lidnum <- fname %>% str_extract("\\-(PHO|LID)\\d\\-") %>% str_replace_all("\\-", "")

  classnum <- fname %>% str_extract("\\-C\\d\\-") %>%  str_replace_all("\\-", "")

  ahdnum <- fname %>% str_extract("AHD_\\d+") %>% str_extract("\\d+") %>% as.numeric()

  # Read spatial extent and acquisition date from html doc
  html <- xml2::read_html(path.html)
  nodes <- rvest::html_nodes(html, "p")

  node.txt <- rvest::html_text(nodes)

  bounds <- node.txt %>% str_subset(ri("bounding"))

  # longitude and latitude - when parsing text, allow for values without decimal point
  lon <- bounds %>% str_subset(ri("longitude")) %>% str_extract("\\-?\\d+(\\.\\d+)?\\s*$") %>% as.numeric()
  lat <- bounds %>% str_subset(ri("latitude")) %>% str_extract("\\-?\\d+(\\.\\d+)?\\s*$") %>% as.numeric()

  epsg <- node.txt %>% str_subset(ri("epsg")) %>% str_extract("\\d+") %>% as.numeric()
  datum <- node.txt %>% str_subset(ri("horizontal datum")) %>% str_extract("[A-Z0-9]+$")
  startdate <- node.txt %>% str_subset(ri("capture start date")) %>% str_extract("[0-9\\-]+") %>% as.Date()
  enddate <- node.txt %>% str_subset(ri("capture end date")) %>% str_extract("[0-9\\-]+") %>% as.Date()

  list(
    filename = fname,
    locality = locality,
    lidnum = lidnum,
    classnum = classnum,
    ahdnum = ahdnum,
    startdate = startdate,
    enddate = enddate,
    bounds.lon = sort(lon),
    bounds.lat = sort(lat),
    epsg = epsg,
    datum = datum)
}


#' Read the header information in a LAS file
#'
#' This function calls \code{\link[rlas]{read.lasheader}} to read header data
#' information from a given LAS file, then returns the data as a nested, named
#' list. The only difference to calling the \code{rlas} function directly is
#' that any problems reading the header will result in a \code{NULL} return
#' value rather than an error, and the element names in the returned list are
#' all lower case with spaces removed to make name matching easier.
#'
#' @param path Path to a LAS file.
#'
#' @param fix.names If TRUE, change all element names in the nested list of
#'   header information to lower case with spaces removed.
#'
#' @return Header information as a nested list
#'
#' @export
#'
read_las_header <- function(path, fix.names = TRUE) {
  if (length(path) != 1) stop("Only a single path is supported at the moment")

  do_fix_names <- function(x) {
    if (is.list(x)) {
      nm <- names(x)
      if (!is.null(nm)) names(x) <- tolower( stringr::str_replace_all(nm, "\\s+", "") )
      for (i in 1:length(x)) {
        x[[i]] <- Recall(x[[i]])
      }
    }
    x
  }

  hdr <- tryCatch(
    rlas::read.lasheader(path),

    error = function(cond) NULL,
    warning = function(cond) NULL
  )

  if (fix.names && !is.null(hdr)) do_fix_names(hdr)
  else hdr
}


.str_length <- function(s) {
  stringr::str_length( stringr::str_trim(s) )
}

.str_ends_with <- function(s, end) {
  stringr::str_detect(s, paste0(end, "$"))
}
