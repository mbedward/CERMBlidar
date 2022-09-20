## template-depr_fun.r: Roxygen template to deprecate functions
#' @name <%= fun %>-deprecated
#' @usage <%= gsub("\n", "\n#' ", roxygen2:::function_usage(fun, formals(fun))) %>
#' @seealso \code{\link{CERMBlidar-deprecated}}
#' @keywords internal
