% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{MCMC}
\alias{MCMC}
\title{Metropolis sampler after a Metropolis within Gibbs Sampler}
\usage{
MCMC(
  niterMwG,
  niterMH,
  parwalk,
  parinit,
  Rexp,
  tdensD,
  alpha,
  parprior,
  adaptive,
  calibration
)
}
\arguments{
\item{niterMwG}{number of iterations of the Metropolis withing Gibbss algorithm mainly used to
infer a covariance matrix for the random walk in the upcoming MH algorithm}

\item{niterMH}{number of iterations of the Metropolis Hastings algorithm}

\item{parwalk}{vector of variances for the random walk}

\item{parinit}{vector of initial parameters rhos, variance of the observation noise,
variance of the discrepancy and (if calibration) thetas}

\item{Rexp}{residual of field data vs simulator}

\item{tdensD}{output of the function tensordist}

\item{alpha}{power in the exponential kernel}

\item{parprior}{prior parameters for the parameters}

\item{adaptive}{boolean to adapt the random walk}

\item{calibration}{list of needed field if a calibration is to be
performed in this order 1) computer model, 2) field observation, 3) field variable, 4) boolean indicating whether calibration
is to be performed or not.}
}
\value{
Posterior sample in the transformed space
}
\description{
Metropolis sampler after a Metropolis within Gibbs Sampler
}
\details{
The MwG sampler is used to estimate a covariance matrix for the random walk in the upcoming Metropolis sampler.
}
