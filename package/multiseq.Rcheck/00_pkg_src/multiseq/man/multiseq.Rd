% Generated by roxygen2 (4.0.2): do not edit by hand
\name{multiseq}
\alias{multiseq}
\title{Estimate underlying signal from count data \code{x} and optionally the effect of a covariate \code{g}.}
\usage{
multiseq(x = NULL, g = NULL, read.depth = NULL, reflect = FALSE,
  baseline = "inter", overall.effect = TRUE, overall.loglr = FALSE,
  cxx = TRUE, maxlogLR = NULL, verbose = FALSE,
  glm.approx.param = list(), ashparam = list(), learn.pi.param = list())
}
\arguments{
\item{x:}{a matrix (or a vector) of \code{nsig} by \code{n} counts where \code{n} should be a power of 2 or a vector of size \code{n}.}

\item{read.depth:}{an \code{nsig}-dimensional vector containing the total number of reads for each sample (used to test for association with the total intensity); defaults to NULL.}

\item{reflect:}{bool, if TRUE signal is reflected, if FALSE signal is not reflected. Defaults to TRUE if n is not power of 2. See \code{\link{reflectSignal}} for details.}

\item{baseline:}{a string, can be "inter" or "grp" or a number. Uses intercept \code{g=0} as baseline ("inter") or the group with the smallest \code{g} as baseline ("grp") or specifies value of \code{g} that should be baseline (number). If center==FALSE and baseline=="inter", then baseline will be overwritten and automatically set to "grp".}

\item{g:}{a \code{nsig}-dimensional vector containing group indicators/covariates for each sample.}

\item{overall.effect:}{bool, indicating whether to include overall mean into effect estimates (TRUE, default) or not (FALSE)}

\item{overall.loglr:}{bool, indicating if multiseq should be used to estimate logLR for overall effect (TRUE) or DESeq2 should be used (FALSE, default).}

\item{cxx:}{bool, indicating whether to use c++ code (faster) (TRUE, default) or R code (FALSE)}

\item{maxlogLR:}{a positive number, defaults to NULL, if \code{maxlogLR} is provided as a positive number, the function returns this number as \code{logLR} when \code{logLR} is infinite.}

\item{verbose:}{bool, defaults to FALSE, if TRUE \code{\link{multiseq}} also outputs \code{logLR$scales} (scales contains (part of) \pkg{ashr} output for each scale), \code{fitted.g}, and \code{fitted.g.intercept}.}

\item{glm.approx.param:}{a list of parameters to be passed to \code{glm.approx}; default values are set by function \code{\link{setGlmApproxParam}}.}

\item{ashparam:}{a list of parameters to be passed to \code{ash}; default values are set by function \code{\link{setAshParam}}.}

\item{learn.pi.param:}{a list of parameters to decide if and how previously learned parameters from \code{ash} should be used; default values are set by function \code{\link{setLearnPiParam}}.}
}
\value{
\code{multiseq} returns an object of \code{\link[base]{class}} "multiseq", a list with the following elements (or a simplified list if \code{verbose=FALSE} or \code{smoothing=FALSE}) \cr
\item{baseline.mean}{an \code{nsig}-vector with the posterior mean of baseline log(intensity)}
\item{baseline.var}{an \code{nsig}-vector with the posterior baseline variance}
\item{effect.mean}{an \code{nsig}-vector with the posterior effect mean}
\item{effect.var}{an \code{nsig}-vector with the posterior effect variance}
\item{logLR}{a list with elements \code{value} specifying the log likelihood ratio, \code{scales} a \code{J+1} vector specifying the logLR at each scale, \code{isfinite} bool specifying if the log likelihood ratio is finite}
\item{fitted.g}{a list of \code{J+1} mixture of normal models fitted using \pkg{ashr}, \code{J=log2(n)}}
\item{fitted.g.intercept}{a list of \code{J} mixture of normal models fitted using \pkg{ashr} on the intercept, \code{J=log2(n)}}
}
\description{
Given a \code{m} by \code{n} matrix of count data \code{X}, we model them as \code{x}_{ib} \sim Poi(\lambda_{ib}) where $log(\lambda_b^i) = \mu_b^o + g^i \beta_b^o + \nu_b^i$ for samples $i=1,...,m$, and $b=1,...,n$, where
$\mu_b^o$ and $\beta_b^o$ are fixed parameters of interest, and $\nu__b^i$ is a random effect modelling extra-Poisson variation. Here the parameters $\mu_b^o$ and
$\beta_b^o$ are assumed to be spatially structured (or smooth if they are treated as being sampled from a continuous function), so the model can be reparametrized as $\alpha_{sl}^i = \mu_{sl} + g^i \beta_{sl} + u_{sl}^i$ using a
multiscale representation, which is the actual model fitted by the function \code{multiseq} (for details see companion paper by Shim et al. (2016)). However, outputs are on the original scale ie \code{multiseq} returns estimates
(posterior means) of $\mu_b^o$ and $\beta_b^o$, along with their posterior standard deviations. Additionally, a log-likelihood ratio for the presence of an "effect" ie deviation of $\beta^o$ from 0 is also included as an output.
}
\examples{
#load data contained in dat
data(dat, package="multiseq")
res <- multiseq(x=dat$x, g=dat$g, minobs=1, lm.approx=FALSE, read.depth=dat$read.depth)
}

