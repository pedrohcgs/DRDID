% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/std_ipw_did_rc.R
\name{std_ipw_did_rc}
\alias{std_ipw_did_rc}
\title{Standardized inverse probability weighted DiD estimator, with repeated cross-section data}
\usage{
std_ipw_did_rc(
  y,
  post,
  D,
  covariates,
  i.weights = NULL,
  boot = FALSE,
  boot.type = "weighted",
  nboot = NULL,
  inffunc = FALSE,
  trim.level = 0.995
)
}
\arguments{
\item{y}{An \eqn{n} x \eqn{1} vector of outcomes from the both pre and post-treatment periods.}

\item{post}{An \eqn{n} x \eqn{1} vector of Post-Treatment dummies (post = 1 if observation belongs to post-treatment period,
and post = 0 if observation belongs to pre-treatment period.)}

\item{D}{An \eqn{n} x \eqn{1} vector of Group indicators (=1 if observation is treated in the post-treatment, =0 otherwise).}

\item{covariates}{An \eqn{n} x \eqn{k} matrix of covariates to be used in the propensity score estimation. Please add a column of ones if you want to include an intercept.
If covariates = NULL, this leads to an unconditional DID estimator.}

\item{i.weights}{An \eqn{n} x \eqn{1} vector of weights to be used. If NULL, then every observation has the same weights. The weights are normalized and therefore enforced to have mean 1 across all observations.}

\item{boot}{Logical argument to whether bootstrap should be used for inference. Default is FALSE.}

\item{boot.type}{Type of bootstrap to be performed (not relevant if \code{boot = FALSE}). Options are "weighted" and "multiplier".
If \code{boot = TRUE}, default is "weighted".}

\item{nboot}{Number of bootstrap repetitions (not relevant if \code{boot = FALSE}). Default is 999.}

\item{inffunc}{Logical argument to whether influence function should be returned. Default is FALSE.}

\item{trim.level}{The level of trimming for the propensity score. Default is 0.995.}
}
\value{
A list containing the following components:
\item{ATT}{The IPW DID point estimate.}
\item{se}{ The IPW DID standard error}
\item{uci}{Estimate of the upper bound of a 95\% CI for the ATT}
\item{lci}{Estimate of the lower bound of a 95\% CI for the ATT}
\item{boots}{All Bootstrap draws of the ATT, in case bootstrap was used to conduct inference. Default is NULL}
\item{att.inf.func}{Estimate of the influence function. Default is NULL}
 \item{call.param}{The matched call.}
 \item{argu}{Some arguments used (explicitly or not) in the call (panel = FALSE, normalized = TRUE, boot, boot.type, nboot, type="ipw")}
}
\description{
\code{std_ipw_did_rc} is used to compute inverse probability weighted (IPW) estimators for the ATT
 in DID setups with stationary repeated cross-sectional data. IPW weights are normalized to sum up to one, that is,
 the estimator is of the Hajek type.
}
\examples{
# use the simulated data provided in the package
covX = as.matrix(cbind(1, sim_rc[,5:8]))
# Implement normalized IPW DID estimator
std_ipw_did_rc(y = sim_rc$y, post = sim_rc$post, D = sim_rc$d,
               covariates= covX)

}
\references{
\cite{Abadie, Alberto (2005), "Semiparametric Difference-in-Differences Estimators",
Review of Economic Studies, vol. 72(1), p. 1-19, \doi{10.1111/0034-6527.00321}.
}


\cite{Sant'Anna, Pedro H. C. and Zhao, Jun. (2020),
"Doubly Robust Difference-in-Differences Estimators." Journal of Econometrics, Vol. 219 (1), pp. 101-122,
\doi{10.1016/j.jeconom.2020.06.003}}
}
