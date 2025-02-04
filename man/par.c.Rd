% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PAR_combine.R
\name{par.c}
\alias{par.c}
\title{Combine parameter blocks}
\usage{
par.c(..., fit = NULL)
}
\arguments{
\item{...}{\code{param_block} arrays and/or character/numeric parameter selections to extract from \code{fit}.
Additional arguments controlling how parameters are extracted from \code{fit} are determined
automatically based on any \code{param_block} arrays in \code{...}. If \code{...} only consist
of parameter selections, the function defaults to extracting posterior samples from \code{fit}. All 
\code{param_block} arrays in \code{...} must have compatible rows and slices (i.e., same number
of chains, same names).
\code{param_block} arrays with a \code{param_type} of "\code{chains}" have rows compatible with any other
\code{param_block} array (they can be coerced to another \code{param_type} on the fly), with the exception
of other \code{chains} \code{param_block} arrays that
have a different number of rows. Otherwise, \code{param_block} arrays generally have compatible rows with other arrays
of the same \code{param_type} provided they share least one row (i.e., quantile or diagnostic) in common.}

\item{fit}{An object of class "\code{evorates_fit}" from which to extract parameters given by the 
character/numeric elements of \code{...}. Ignored if \code{NULL} (the default).}
}
\value{
An array of class "\code{param_block}" with a \code{param_type} determined by the elements
of \code{...} (tends to default "\code{chains}" when in doubt). The dimensions of these arrays will generally
go in the order of iterations/quantiles/diagnostics, then parameters, then chains.  Any dimensions of length 1 are
collapsed and stored as attributes. The resulting array contains all parameters specified by \code{...} as columns,
ordered as they are in \code{...}.
}
\description{
This function is a wrapper for the \code{c()} method for \code{param_block} arrays that allows you
to also specify a fitted evorates model from which to extract parameters specified via \code{...}.
}
\details{
I had to create a separate function for this because R does not default to the \code{c()} method for
\code{param_block} arrays if character or numeric arguments are passed to \code{c()}. This function
simply ensures R calls the proper method no matter what's passed to \code{...}.
}
\examples{
#get whale/dolphin evorates fit
data("cet_fit")

#how are average rates for some focal clades affected by R_sig2/R_mu estimates?
parblock <- get.bg.rate(fit = cet_fit,
                        node.groups = setNames(list('Mesoplodon','Orcinus',c('Pseudorca','Feresa')),
                                               c('Mesoplodon','Orca','Globicephalinae')),
                        )
parblock <- par.c("R_mu", "R_sig2", parblock, fit = cet_fit)
plot(parblock \%chains\% "Mesoplodon" ~ parblock \%chains\% "R_sig2"
plot(parblock \%chains\% "Orca" ~ parblock \%chains\% "R_mu")

#automatic conversion based on param_block type
parblock <- get.bg.rate(fit = cet_fit,
                        node.groups = setNames(list('Mesoplodon','Orcinus',c('Pseudorca','Feresa')),
                                               c('Mesoplodon','Orca','Globicephalinae')),
                        type = "quantiles")
par.c("R_mu", "R_sig2", parblock, fit = cet_fit)

#can mix numeric and character selections
par.c("R_0", "R_mu", 1, 28, fit = cet_fit)


}
\seealso{
\link{param_block-class} for general information on \code{param_block} arrays and 
\code{\link[=grapes-chains-grapes]{\%chains\%}()},
\code{\link[=grapes-quantiles-grapes]{\%quantiles\%}()},
\code{\link[=grapes-means-grapes]{\%means\%}()},
\code{\link[=grapes-diagnostics-grapes]{\%diagnostics\%}()},
and \code{\link[=grapes-select-grapes]{\%select\%}()} for more information on
\code{param_block} operators.
}
