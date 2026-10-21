context("Influence function of the DR DiD estimators with repeated cross sections")

# att.inf.func should be the derivative of the ATT with respect to the weight of
# each observation. Giving observation i weight 1 + h (all others weight 1) moves
# the empirical distribution by eps(h) = h / (n + h) towards a point mass at i, so
#   (ATT(h) - ATT(-h)) / (eps(h) - eps(-h))  ->  att.inf.func[i],
# with an error of order h^2 (about 1e-8 relative here for h = 0.01).
# The check is done on observations from every (D, post) cell, including the
# pre-period comparison units, whose influence values carry the estimation effect
# of the pre-period outcome regression.

make_rc_inffunc_data <- function() {
  set.seed(20261021)
  n <- 400
  x <- stats::rnorm(n)
  # the logit in x is misspecified, so the weighted covariate means of the
  # treated and comparison cells differ
  D <- as.numeric(stats::runif(n) < stats::plogis(-0.5 + 1.5 * pmax(x, 0)^2))
  post <- as.numeric(stats::runif(n) < 0.5)
  y <- 1 + 2 * x + 0.5 * D + post * (0.5 + x + D) +
    stats::rnorm(n) * ifelse(D == 0 & post == 0, 3, 1)
  list(y = y, post = post, D = D, X = cbind(1, x), n = n)
}

numerical_inffunc <- function(estimator, d, i, h = 0.01) {
  att_w <- function(hh) {
    w <- rep(1, d$n)
    w[i] <- 1 + hh
    estimator(d$y, d$post, d$D, d$X, i.weights = w)$ATT
  }
  eps <- function(hh) hh / (d$n + hh)
  (att_w(h) - att_w(-h)) / (eps(h) - eps(-h))
}

test_that("drdid_rc and drdid_rc1 influence functions match numerical derivatives", {
  d <- make_rc_inffunc_data()
  # two observations from each (D, post) cell
  cells <- interaction(d$D, d$post)
  idx <- unname(unlist(lapply(split(seq_len(d$n), cells), utils::head, 2)))

  for (estimator in list(drdid_rc, drdid_rc1)) {
    inf_func <- estimator(d$y, d$post, d$D, d$X, inffunc = TRUE)$att.inf.func
    num_inf_func <- vapply(idx, function(i) numerical_inffunc(estimator, d, i),
                           numeric(1))
    expect_equal(as.vector(inf_func[idx]), num_inf_func, tolerance = 1e-6)
  }
})
