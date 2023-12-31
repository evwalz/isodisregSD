% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/partialOrders.R
\name{compOrd_cpp}
\alias{compOrd_cpp}
\title{SD order relation}
\usage{
compOrd_cpp(X)
}
\arguments{
\item{X}{a numeric matrix or a data frame containing numeric or ordered
factor variables with at least two columns and for which there are no distinct
rows that are exactly equal.}
}
\value{
A list containing

\item{\code{paths}}{a two-column matrix giving all pairs of indices
\code{(i,j)} which satisfy \code{(X[i, ] <= X[j, ])}.}

\item{\code{colOrder}}{a matrix of the columnwise orders of \code{X}. Used to
compute paths and required for other function calls.}
}
\description{
Compares all rows of a numeric matrix or data frame \code{X}
with respect to the SD order.
}
\keyword{internal}
