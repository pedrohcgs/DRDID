# Validate before mean normalization can conceal negative input weights.
validate_sampling_weights <- function(i.weights, n) {
  if (!(is.numeric(i.weights) || is.logical(i.weights)) ||
      length(i.weights) != n || any(!is.finite(i.weights)) ||
      any(i.weights < 0) || !any(i.weights > 0) ||
      !is.finite(mean(i.weights)) || mean(i.weights) <= 0) {
    stop("i.weights must be finite, non-negative, have positive mean, and contain one value per observation.",
         call. = FALSE)
  }
  invisible(NULL)
}
