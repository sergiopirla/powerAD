% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/powerADgen.R
\name{powerADgen}
\alias{powerADgen}
\title{powerADgen}
\usage{
powerADgen(
  individuals,
  samples,
  metric,
  r,
  p.value = 0.05,
  min,
  max,
  data,
  id,
  affect,
  time,
  perm = 100
)
}
\arguments{
\item{individuals}{Number of individuals sampled.}

\item{samples}{Number of affect observations per individual.}

\item{metric}{Metric of interest. Must be a single character from "Average", "Rel.SD", "SD", "RMSSD", "TKEO", "PAC" or "Autocorrelation".}

\item{r}{Number from 0.01 to 0.99 indicating the expected effect size (Pearson correlation) of interest.}

\item{p.value}{Alpha level. Must be one of the following numbers: 0.01, 0.05, 0.001, 0.005, 0.001. Defaults to 0.05.}

\item{min}{Affect scale lower bound.}

\item{max}{Affect scale upper bound.}

\item{data}{Data frame containing affect observations used perform power analysis.}

\item{id}{Name of variable containing id observations.}

\item{affect}{Name of variable containing affect observations.}

\item{time}{Name of variable containing observation order/time stamp of affect reports.}

\item{perm}{Number of simulations employed in power analysis. For robust power calculations, we recommend using a minimum of 1000 simulations.}
}
\description{
Estimates the statistical power of an affect dynamics time series study for a given sampling, effect size, alpha level and metric of interest.
}
\references{
Pirla, Taquet and Quoidbach (2021). ADD REFERENCE
}
\keyword{Affect}
\keyword{Analysis}
\keyword{Dynamics,}
\keyword{Power}
