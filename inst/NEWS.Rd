\name{NEWS}
\title{News for Package 'systemicrisk'}

\section{Changes in version 0.4.2}{
  \itemize{
    \item Bug fix for matrices of dimension > 1024.
    \item Added specific support for some non-square matrices: Model.lambda.constant.nonsquare,
  Model.p.constant.nonsquare, calibrate_ER.nonsquare.
  }
}


\section{Changes in version 0.4.1}{
  \itemize{
    \item Addjusted that non-square matrices can be used (using fixed
  matrices for p and lambda). See the new vignette as an example.
  }
}


\section{Changes in version 0.4}{
  \itemize{
    \item Added methods to calibrate models to a given density.
  }
}

\section{Changes in version 0.3.1}{
  \itemize{
    \item Fixed parameter range checks for Fitness models.
  }
}
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