# Obsolete or unused functions

# Metrics function for summary statistics
#
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


# Calculate highest density interval for a vector of values
#
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

