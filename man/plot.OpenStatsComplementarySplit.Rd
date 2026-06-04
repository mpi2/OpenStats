\name{plot.OpenStatsComplementarySplit}
\alias{plot.OpenStatsComplementarySplit}
\title{plot for an `OpenStatsComplementarySplit` object}
\description{This function visualises an `OpenStatsComplementarySplit` object}

\usage{
 \method{plot}{OpenStatsComplementarySplit}(x, main = "Final Model", ask = FALSE, mfrow = c(2, 2), ...)
}
\arguments{
  \item{x}{
  an instance of `OpenStatsComplementarySplit` result 
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
The plot function creates some visualisations for the split results from the linear mixed model framework. Each level of partitioning variables (see `variables` in the `OpenStatsComplementarySplit` function manual) produces a set of plots listed below: \cr

- Residual versus fitted values \cr 
- Residual density plot and the normality test p-value\cr
- Residual Q-Q plot\cr
- The density plot of the response variable and the normality test p-value\cr
}

\value{
Not applicable
}

\seealso{
\code{\link{OpenStatsComplementarySplit}}, \code{\link{print.OpenStatsComplementarySplit}}, \code{\link{summary.OpenStatsComplementarySplit}}
}

\examples{
example(OpenStatsComplementarySplit)
}

