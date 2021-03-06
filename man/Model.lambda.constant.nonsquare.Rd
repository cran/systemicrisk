% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HierarchicalModels.R
\name{Model.lambda.constant.nonsquare}
\alias{Model.lambda.constant.nonsquare}
\title{Model for a Constant lambda and Non-Square Matrices}
\usage{
Model.lambda.constant.nonsquare(lambda, nrow, ncol)
}
\arguments{
\item{lambda}{paramer for the size of the liabilities. A single numeric value.}

\item{nrow}{number of rows of the matrix.}

\item{ncol}{number of columns of the matrix.}
}
\value{
the resulting model.
}
\description{
This model assumes that the parameter lambda is known.
}
\examples{
m <- Model.lambda.constant.nonsquare(nrow=5,ncol=3,lambda=0.25)
m$matr(m$rtheta())
}
