% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/psis.R, R/sis.R, R/tis.R
\name{psis}
\alias{psis}
\alias{psis.array}
\alias{psis.matrix}
\alias{psis.default}
\alias{is.psis}
\alias{is.sis}
\alias{is.tis}
\title{Pareto smoothed importance sampling (PSIS)}
\usage{
psis(log_ratios, ...)

\method{psis}{array}(log_ratios, ..., r_eff = 1, cores = getOption("mc.cores", 1))

\method{psis}{matrix}(log_ratios, ..., r_eff = 1, cores = getOption("mc.cores", 1))

\method{psis}{default}(log_ratios, ..., r_eff = 1)

is.psis(x)

is.sis(x)

is.tis(x)
}
\arguments{
\item{log_ratios}{An array, matrix, or vector of importance ratios on the log
scale (for PSIS-LOO these are \emph{negative} log-likelihood values). See the
\strong{Methods (by class)} section below for a detailed description of how
to specify the inputs for each method.}

\item{...}{Arguments passed on to the various methods.}

\item{r_eff}{Vector of relative effective sample size estimates containing
one element per observation. The values provided should be the relative
effective sample sizes of \code{1/exp(log_ratios)} (i.e., \code{1/ratios}).
This is related to the relative efficiency of estimating the normalizing
term in self-normalizing importance sampling. If \code{r_eff} is not
provided then the reported PSIS effective sample sizes and Monte Carlo
error estimates can be over-optimistic. If the posterior draws are (near)
independent then \code{r_eff=1} can be used. \code{r_eff} has to be a scalar (same
value is used for all observations) or a vector with length equal to the
number of observations. The default value is 1. See the \code{\link[=relative_eff]{relative_eff()}}
helper function for computing \code{r_eff}.}

\item{cores}{The number of cores to use for parallelization. This defaults to
the option \code{mc.cores} which can be set for an entire R session by
\code{options(mc.cores = NUMBER)}. The old option \code{loo.cores} is now
deprecated but will be given precedence over \code{mc.cores} until
\code{loo.cores} is removed in a future release. \strong{As of version
2.0.0 the default is now 1 core if \code{mc.cores} is not set}, but we
recommend using as many (or close to as many) cores as possible.
\itemize{
\item Note for Windows 10 users: it is \strong{strongly}
\href{https://github.com/stan-dev/loo/issues/94}{recommended} to avoid using
the \code{.Rprofile} file to set \code{mc.cores} (using the \code{cores} argument or
setting \code{mc.cores} interactively or in a script is fine).
}}

\item{x}{For \code{is.psis()}, an object to check.}
}
\value{
The \code{psis()} methods return an object of class \code{"psis"},
which is a named list with the following components:

\describe{
\item{\code{log_weights}}{
Vector or matrix of smoothed (and truncated) but \emph{unnormalized} log
weights. To get normalized weights use the
\code{\link[=weights.importance_sampling]{weights()}} method provided for objects of
class \code{"psis"}.
}
\item{\code{diagnostics}}{
A named list containing two vectors:
\itemize{
\item \code{pareto_k}: Estimates of the shape parameter \eqn{k} of the
generalized Pareto distribution. See the \link{pareto-k-diagnostic}
page for details.
\item \code{n_eff}: PSIS effective sample size estimates.
}
}
}

Objects of class \code{"psis"} also have the following \link[=attributes]{attributes}:
\describe{
\item{\code{norm_const_log}}{
Vector of precomputed values of \code{colLogSumExps(log_weights)} that are
used internally by the \code{weights} method to normalize the log weights.
}
\item{\code{tail_len}}{
Vector of tail lengths used for fitting the generalized Pareto distribution.
}
\item{\code{r_eff}}{
If specified, the user's \code{r_eff} argument.
}
\item{\code{dims}}{
Integer vector of length 2 containing \code{S} (posterior sample size)
and \code{N} (number of observations).
}
\item{\code{method}}{
Method used for importance sampling, here \code{psis}.
}
}
}
\description{
Implementation of Pareto smoothed importance sampling (PSIS), a method for
stabilizing importance ratios. The version of PSIS implemented here
corresponds to the algorithm presented in Vehtari, Simpson, Gelman, Yao,
and Gabry (2024).
For PSIS diagnostics see the \link{pareto-k-diagnostic} page.
}
\section{Methods (by class)}{
\itemize{
\item \code{psis(array)}: An \eqn{I} by \eqn{C} by \eqn{N} array, where \eqn{I}
is the number of MCMC iterations per chain, \eqn{C} is the number of
chains, and \eqn{N} is the number of data points.

\item \code{psis(matrix)}: An \eqn{S} by \eqn{N} matrix, where \eqn{S} is the size
of the posterior sample (with all chains merged) and \eqn{N} is the number
of data points.

\item \code{psis(default)}: A vector of length \eqn{S} (posterior sample size).

}}
\examples{
log_ratios <- -1 * example_loglik_array()
r_eff <- relative_eff(exp(-log_ratios))
psis_result <- psis(log_ratios, r_eff = r_eff)
str(psis_result)
plot(psis_result)

# extract smoothed weights
lw <- weights(psis_result) # default args are log=TRUE, normalize=TRUE
ulw <- weights(psis_result, normalize=FALSE) # unnormalized log-weights

w <- weights(psis_result, log=FALSE) # normalized weights (not log-weights)
uw <- weights(psis_result, log=FALSE, normalize = FALSE) # unnormalized weights



}
\references{
Vehtari, A., Gelman, A., and Gabry, J. (2017). Practical Bayesian model
evaluation using leave-one-out cross-validation and WAIC.
\emph{Statistics and Computing}. 27(5), 1413--1432. doi:10.1007/s11222-016-9696-4
(\href{https://link.springer.com/article/10.1007/s11222-016-9696-4}{journal version},
\href{https://arxiv.org/abs/1507.04544}{preprint arXiv:1507.04544}).

Vehtari, A., Simpson, D., Gelman, A., Yao, Y., and Gabry, J. (2024).
Pareto smoothed importance sampling. \emph{Journal of Machine Learning Research},
25(72):1-58.
\href{https://jmlr.org/papers/v25/19-556.html}{PDF}
}
\seealso{
\itemize{
\item \code{\link[=loo]{loo()}} for approximate LOO-CV using PSIS.
\item \link{pareto-k-diagnostic} for PSIS diagnostics.
\item The \strong{loo} package \href{https://mc-stan.org/loo/articles/index.html}{vignettes}
for demonstrations.
\item The \href{https://mc-stan.org/loo/articles/online-only/faq.html}{FAQ page} on
the \strong{loo} website for answers to frequently asked questions.
}
}
