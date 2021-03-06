% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/crossover_func.R
\name{crossover_func}
\alias{crossover_func}
\title{Crossover Function}
\usage{
crossover_func(crossover_matrix_row, Kvar, Kint, pc)
}
\arguments{
\item{crossover_matrix_row}{Vector containing both parent models.}

\item{Kvar}{Maximum number of variables allowed in the final model}

\item{Kint}{Maximum number of interactions allowed in the final model}

\item{pc}{Legacy variable (equal to 0.5)}
}
\value{
Vector containing two children models.
}
\description{
Performs crossover between two models. Internal use only.
}
\examples{
p1 <- c("x1","0","x4","x3","0","x8","1","0","0","0","0","1")
p2 <- c("x1","x2","x4","x7","0","x8","1","0","1","1","0","0")
crossover_matrix_row <- c(p1,p2)
Kvar <- 6
Kint <- 0
pc <- 0.5

\dontrun{
crossover_func(crossover_matrix_row,Kvar,Kint,pc)
}

}
\keyword{internal}
