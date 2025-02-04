% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/OUT_OPS_main.R
\name{\%quantiles\%}
\alias{\%quantiles\%}
\alias{\%q\%}
\title{Extract posterior distribution quantiles from a fitted evorates model}
\usage{
fit \%quantiles\% select

fit \%q\% select
}
\arguments{
\item{fit}{An object of class "\code{evorates_fit}" or "\code{param_block}". If \code{fit} is a
\code{param_block} array, it must have a \code{param_type} of "\code{chains}" or
"\code{quantiles}".}

\item{select}{A list with two elements (2nd element is optional):
\itemize{
\item{A character or numeric vector for selecting parameters. If a character vector, entries are
matched to parameter names using regular expressions, with one key exception: if any entries are
\emph{exact} matches to sampling parameter names, these will be used to select parameters from
the \code{sampler.params} element of \code{fit} instead, if it exists (see details). If a 
numeric vector, entries are matched to edge indices to select branchwise rate parameters;
these can be negative to instead \emph{exclude} branchwise rates. If no branchwise rate parameters
are found, a number \code{i} will instead select the \code{i}th parameter in \code{fit}.}
\item{A numeric, integer, or character vector for specifying quantiles. If unsupplied, default quantiles
are selected. Default quantiles are are set to the quantiles already stored in \code{fit}, if they exist,
otherwise they are set to the lower bound of the 95\% credible interval, median, and upper bound (i.e., 2.5\%,
50\%, and 97.5\%). If a numeric vector, this will simply be the quantile to extract (i.e., 0.1 would be
the 10\% quantile). If an integer vector, an integer \code{i} will extract the \code{i}th default quantile.
If a character vector, numbers as taken as percents (i.e., 10\% or even 10 would correspond to 0.1).
Quantiles greater than 100\% or less than 0\% are rounded to 100\% and 0\%, respectively.}
}}
}
\value{
An array of class "\code{param_block}" with a \code{param_type} of "\code{quantiles}".
The dimension of these arrays will generally go in the order of quantiles, then parameters,
then chains. Any dimensions of length 1 are collapsed and stored as attributes. If \code{fit}
is a \code{chains param_block} array, parameters are automatically renamed
"\code{quantiles(<parameter name>)}" to help keep track of parameter manipulations.
}
\description{
This operator extracts quantiles from the posterior distributions for particular parameters from
an \code{evorates_fit} object or \code{param_block} array. This can be used to, for example, get
posterior medians or credible intervals.
}
\details{
In the case that a numeric vector is provided to select parameters and no parameters
with names following the pattern "\code{R_i}" are found, the function then looks for the pattern
"\code{Rdev_i}", then "\code{uncent_Rdev_i}". If neither of these are found, then it finally
defaults to selecting the \code{i}th parameters. If a single parameter name involves multiple
"\code{R_i}" patterns, the first pattern is always used (e.g., \code{R_1\%-\%R_2} would correspond
to \code{1}). Similarly, if multiple parameters match to a single number, the first parameter is
always used (similar to behavior of \link{match}).

The \code{sampler.params} element of \code{fit} always includes 9 parameters, named:
"\code{accept_stat__}", "\code{treedepth__}", "\code{stepsize__}", "\code{divergent__}",
"\code{n_leapfrog__}", "\code{energy__}", "\code{prior__}", "\code{lik__}", and "\code{post__}".
This tends to include both warmup and non-warmup samples, but warmup samples are
automatically excluded if both sampler parameters and normal parameters are selected by this
function. Most of these parameters are rather technical quantities used to tune the
Hamiltonian Monte Carlo sampler run by Stan (see Stan manual for further details on what they
mean). Generally, users will only want to look at the last 3 parameters, which give the
(log) prior probability, likelihood, and posterior probability, respectively, of sampled
parameters. Note that these are on the sampling scale and will differ from those for the
originally-scaled data by a constant. Also, the posterior probability will be affected
by what \code{lik.power} was set to for fitting the model. In some cases, users may also
wish to look at what parameter values are associated with divergent transitions (i.e.,
iterations where \code{divergent__ = 1}), which indicate regions of parameter space where
the sampler got "stuck", yielding potentially misleading posterior distribution estimates.

This function uses a custom quantiles estimation function that yields the same results as
\link{quantile} with \code{type = 7}. \code{NA}s are always ignored when estimating
quantiles, unless all samples for a parameter are \code{NA}, in which case all
quantile estimates will also be \code{NA}.
}
\examples{
#get whale/dolphin evorates fit
data("cet_fit")

#extracting directly from evorates fit
cet_fit \%quantiles\% "R_mu"
#regular expressions
cet_fit \%quantiles\% "R"
#using . is a quick way to extract ALL parameters!
cet_fit \%quantiles\% "."
#numeric index-based selection
cet_fit \%quantiles\% 1
cet_fit \%quantiles\% -1
#select particular quantiles
#e.g., medians
cet_fit \%quantiles\% list("R_mu", 0.5)
cet_fit \%quantiles\% list("R_mu", 2L)
cet_fit \%quantiles\% list("R_mu", "50\%")
#e.g., 95\% credible intervals
cet_fit \%quantiles\% list("R_mu", c(0.025, 0.975))
cet_fit \%quantiles\% list("R_mu", c(1L, 3L))
cet_fit \%quantiles\% list("R_mu", c("2.5\%", "97.5\%"))
#getting sampler parameters
cet_fit \%quantiles\% "lik__"
#note warmup samples automatically excluded from "lik__" if combined with "R_mu"
cet_fit \%quantiles\% c("R_mu", "lik__")

#extracting from a param_block array
par <- get.bg.rate(fit = cet_fit,
                   node.groups = setNames(list('Mesoplodon','Orcinus',c('Pseudorca','Feresa')),
                                          c('Mesoplodon','Orca','Globicephalinae')),
                   )
par \%quantiles\% list("Mesoplodon", 0.5)
#note change in numeric index behavior
par \%quantiles\% 1


}
\seealso{
Other param_block operators: 
\code{\link{\%chains\%}()},
\code{\link{\%diagnostics\%}()},
\code{\link{\%means\%}()},
\code{\link{\%select\%}()}
}
\concept{param_block operators}
