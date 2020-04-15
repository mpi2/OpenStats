\name{plot.OpenStatsMM}
\alias{plot.OpenStatsMM}
\title{plot for an `OpenStatsMM` object}
\description{This function visualises an `OpenStatsMM` object}

\usage{
 \method{plot}{OpenStatsMM}(x, main = "Final Model", ask = FALSE, mfrow = c(2, 2), ...)
}
\arguments{
  \item{x}{
  an instance of `OpenStatsMM` result from OpenStatsAnalysis(method = 'MM') function
}
  \item{main}{
  a string to be pasted to the title of the plots
}
  \item{ask}{
  see `ask` in `par()` function. Default FALSE
}
  \item{mfrow}{
  the screen partition. see `mfrow` argument in the `par` function. Default c(2,2) then all plots display in one screen.
}
  \item{\dots}{
  other parameters that can be passed to the `plot` function
}
}
\details{
The plot function creates some visualisations for the linear mixed model framework: \cr

- Residual versus fitted values \cr 
- Residual density plot and the normality test p-value\cr
- Residual Q-Q plot\cr
- The density plot of the response variable and the normality test p-value\cr
}

\author{
  Hamed Haseli Mashhadi <hamedhm@ebi.ac.uk>
}

\seealso{
\code{\link{OpenStatsAnalysis}}, \code{\link{plot.OpenStatsFE}}, \code{\link{plot.OpenStatsRR}}
}

\examples{
example(OpenStatsAnalysis)
}

