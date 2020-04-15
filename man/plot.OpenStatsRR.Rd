\name{plot.OpenStatsRR}
\alias{plot.OpenStatsRR}
\title{plot for an `OpenStatsRR` object}
\description{This function visualises an `OpenStatsRR` object}

\usage{
 \method{plot}{OpenStatsRR}(x, main = "Mosaic plot", ask = FALSE, mfrow = c(2, 2), ...)
}
\arguments{
  \item{x}{
  an instance of `OpenStatsRR` result from OpenStatsAnalysis(method = 'RR') function
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
The plot function creates some visualisations for the reference range plus framework \cr

- Mosaic plot of the discretised response versus Genotype/Sex/LifeStage (if they exist in the data\cr
- Mosaic plot of the Sex versus Genotype (if they exist in the data)\cr
}

\author{
  Hamed Haseli Mashhadi <hamedhm@ebi.ac.uk>
}

\seealso{
\code{\link{OpenStatsAnalysis}}, \code{\link{plot.OpenStatsFE}}, \code{\link{plot.OpenStatsMM}}
}

\examples{
example(OpenStatsAnalysis)
}

