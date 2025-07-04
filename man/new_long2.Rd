% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dataset.R
\docType{data}
\name{new_long2}
\alias{new_long2}
\title{longitudinal- survival dataset}
\format{
a tibble of 13 columns and 5639 observations,
\describe{
\item{id}{id value for subjects}
\item{status}{survival status}
\item{time}{ survival time}
\item{y}{longitudinal marker}
\item{visit}{visit time of longitudinal measurements}
\item{x1,x2,...,x7}{ different numeric and categorical variable}

}
}
\usage{
data(new_long2)
}
\description{
A longitudinal dataset with single marker , with different numeric and categorical covariate
}
\keyword{datasets}
