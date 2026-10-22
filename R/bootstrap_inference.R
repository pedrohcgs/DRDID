# Summarize centered bootstrap errors without changing positive-scale inference.
bootstrap_inference <- function(errors, estimate) {
  if (length(estimate) != 1L || !is.finite(estimate) ||
      !is.numeric(errors) || !length(errors) || any(!is.finite(errors))) {
    stop("Bootstrap inference requires a finite estimate and finite bootstrap errors.",
         call. = FALSE)
  }

  se <- stats::IQR(errors) / (stats::qnorm(0.75) - stats::qnorm(0.25))
  if (!is.finite(se)) {
    stop("Bootstrap inference has a non-finite IQR scale.", call. = FALSE)
  }
  if (se == 0) {
    if (all(errors == 0)) {
      return(list(se = 0, lci = estimate, uci = estimate))
    }
    stop(paste0("Bootstrap errors have zero IQR but are not all zero; ",
                "IQR-based bootstrap inference is undefined."), call. = FALSE)
  }

  # Keep the original arithmetic (including quantile names) for positive scales.
  cv <- stats::quantile(abs(errors / se), probs = 0.95)
  lci <- estimate - cv * se
  uci <- estimate + cv * se
  if (!is.finite(cv) || !is.finite(lci) || !is.finite(uci)) {
    stop("Bootstrap inference produced a non-finite critical value or interval.",
         call. = FALSE)
  }
  list(se = se, lci = lci, uci = uci)
}
