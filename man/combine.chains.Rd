% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/OUT_chain.management.R
\name{combine.chains}
\alias{combine.chains}
\title{Combine chains in a fitted evorates model}
\usage{
combine.chains(fit, simplify = TRUE)
}
\arguments{
\item{fit}{An object of class "\code{evorates_fit}".}

\item{simplify}{\code{TRUE} or \code{FALSE}: should the resulting \code{param_block} array(s) in
\code{fit} be  simplified? If \code{TRUE} (the default), dimensions of length 1 in the result are
automatically collapsed, with corresponding information stored as attributes (this is the default
behavior of param_block operators).}
}
\value{
An object of class "\code{evorates_fit}". All previously-existing \code{param_block} arrays stored
in \code{fit} will be included.
}
\description{
This function combines all chains in fitted evorates model sequentially into a single, larger chain.
Generally, this function should only be run after confirming that (via Rhat diagnostics) that each
chain adequately converged!
}
\details{
Chains are combined sequentially, meaning that the beginning of the second chain will follow
the end of the first chain, the third will follow the second, and so on--no permutation business is
attempted. The resulting chain name will simply be all the previous chain names separated by commas
(e.g., "\code{chain 1, chain 2, chain 3, <...>}"). Because of this, initial values in the parameter
diagnostics param_block (if present) are set to \code{NA} and any remaining "warmup iterations"
(defined as iterations present in \code{fit$sampler.params} but not in \code{fit$chains}) are
removed, with all retained iterations reclassified (perhaps misleadingly in some cases)
as non-warmup iterations.
}
\examples{
#get whale/dolphin evorates fit
data("cet_fit")

#make sure all chains co


}
\seealso{
Other chain management: 
\code{\link{exclude.warmup}()},
\code{\link{select.chains}()},
\code{\link{thin.chains}()}
}
\concept{chain management}
