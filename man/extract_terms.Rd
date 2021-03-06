% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract_terms.R
\name{extract_terms}
\alias{extract_terms}
\title{Extract Terms from Model}
\usage{
extract_terms(individual, Kvar, Kint)
}
\arguments{
\item{individual}{Vector containing a model}

\item{Kvar}{Maximum number of variables allowed in the final model}

\item{Kint}{Maximum number of interactions allowed in the final model}
}
\value{
A list specifying which variables are included linearly and which are included nonparametrically.
}
\description{
Extracts linear and nonparametric variables from a model. Internal use only.
}
\keyword{internal}
