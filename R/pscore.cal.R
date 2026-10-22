###################################################################################
# Compute propensity scores using IPT. The logit fit supplies starting values only.
# A returned flag identifies the validated candidate: 0 = trust, 1 = secondary IPT.

pscore.cal <- function(D, int.cov, i.weights, n, bootstrap.index = NULL){
  fail <- function(){
    message <- paste0("IPT propensity-score fit did not converge. Consider explicitly using ",
                      "the traditional DR estimator by rerunning `drdid(..., estMethod = \"trad\")`.")
    if (!is.null(bootstrap.index)) {
      message <- paste0(message, "\nWeighted-bootstrap replicate ", bootstrap.index, " failed.")
    }
    stop(message, call. = FALSE)
  }

  if (any(!is.finite(i.weights)) || any(i.weights < 0) || !any(i.weights > 0)) fail()
  # Common weight scale cannot change an IPT solution, but can change optimizer
  # stopping criteria. Use mean-one weights for both solvers and the initializer.
  fitting.weights <- i.weights / mean(i.weights)
  if (any(!is.finite(fitting.weights))) {
    fitting.weights <- i.weights / max(i.weights)
    fitting.weights <- fitting.weights / mean(fitting.weights)
  }

  # A failed initializer must not substitute an MLE for an IPT solution.
  pslogit <- tryCatch(suppressWarnings(fastglm_fit(int.cov, D, stats::binomial(),
                                                  fitting.weights, method = 3)),
                      error = function(e) NULL)
  init.gamma <- pslogit$coefficients
  if (length(init.gamma) != ncol(int.cov) || any(!is.finite(init.gamma))) {
    init.gamma <- rep(0, ncol(int.cov))
  }

  # The optimizer changes only gamma. Cache the filtered problem and immutable
  # range facts once; direct calls to either criterion prepare their own state.
  prepared <- ipt.loss.prepare(D, int.cov, fitting.weights)
  primary <- tryCatch(suppressWarnings(trust::trust(loss.ps.cal, parinit = init.gamma,
                        rinit = 1, rmax = 1000, iterlim = 1000,
                        D = D, int.cov = int.cov, iw = fitting.weights, .ipt = prepared)),
                      error = function(e) NULL)
  gamma.cal <- primary$argument
  flag <- 0

  # Optimizer status alone does not establish original exponential-IPT balance.
  if (!ipt.valid(gamma.cal, D, int.cov, i.weights)) {
    secondary <- tryCatch(suppressWarnings(stats::nlm(loss.ps.IPT, init.gamma,
                          iterlim = 10000, gradtol = 1e-08, check.analyticals = FALSE,
                          D = D, int.cov = int.cov, iw = fitting.weights, n = n,
                          .ipt = prepared)),
                        error = function(e) NULL)
    gamma.cal <- secondary$estimate
    if (!ipt.valid(gamma.cal, D, int.cov, i.weights)) fail()
    flag <- 1
    warning(paste0("trust algorithm did not converge to a verified IPT fit. ",
                   "Used the secondary IPT algorithm."), call. = FALSE)
  }

  pscore <- as.numeric(stats::plogis(int.cov %*% gamma.cal))
  if (any(!is.finite(pscore))) fail()
  list(pscore = pscore, flag = flag)
}
