\name{NEWS}
\title{News for Package 'systemicrisk'}

\section{Changes in version 0.3 (2015-12-21)}{
  \itemize{
    \item Added Hierarchical Models.
    \item Arbitrary values of the Matrix can now be force to equal 0
  (not just the diagonal). An error message will be produced if no
  matrix exists for the given row and column sums. The function
  findFeasibleMatrix() supersedes getfeasibleMatrix()
    \item Added automatic choice of thinning to achieve desired relative
    effective sample size (choosethin)
    \item Added possibility to fix elements of L to some values.
    \item Added a fitness based model for the matrix that incorporates
  power laws.
  }
}
\section{Changes in version 0.2 (2015-03-20)}{
  \itemize{
    \item First CRAN upload.
  }
}