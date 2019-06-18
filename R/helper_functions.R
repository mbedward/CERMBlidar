# Check that an argument is a scalar. If it is a vector, issue
# a warning and return just the first value
.as_scalar <- function(x) {
  if (length(x) == 1) x
  else {
    nm <- deparse(substitute(x))
    if (length(x) == 0) {
      stop("Expected a value for ", nm)
    }
    else {
      warning("Expected a single value for ", nm, ". Only using first value.")
      x[1]
    }
  }
}
