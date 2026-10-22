# NULL retains the low-level API default; inactive bootstrap never calls this.
validate_bootstrap_count <- function(nboot) {
  if (is.null(nboot)) return(999)
  if (!is.numeric(nboot) || is.complex(nboot) || length(nboot) != 1L || !is.finite(nboot) ||
      nboot < 2 || nboot != floor(nboot)) {
    stop("nboot must be a finite scalar integer of at least 2 when boot = TRUE.",
         call. = FALSE)
  }
  nboot
}
