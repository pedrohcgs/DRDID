# Missing wrapper rows are handled by preprocessing; low-level inputs are complete.
validate_treatment_indicator <- function(D, allow_missing = FALSE) {
  if (!(is.numeric(D) || is.logical(D)) ||
      (!allow_missing && anyNA(D)) ||
      any(!is.na(D) & !(D %in% c(0, 1)))) {
    stop("D must contain only 0 (comparison) and 1 (treated).", call. = FALSE)
  }
  invisible(NULL)
}
