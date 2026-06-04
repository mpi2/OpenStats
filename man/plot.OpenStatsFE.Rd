\name{plot.OpenStatsFE}
\alias{plot.OpenStatsFE}
\title{plot for an `OpenStatsFE` object}
\description{This function visualises an `OpenStatsFE` object}

\usage{
 \method{plot}{OpenStatsFE}(x, main = "Mosaic plot", ask = FALSE, mfrow = c(2, 2), ...)
}
\arguments{
  \item{x}{
  an instance of `OpenStatsFE` result from OpenStatsAnalysis(method = 'FE') function
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

\value{
Not applicable
}

\details{
The plot function creates some visualisations for the Fisher's exact test framework: \cr

- Mosaic plot of the response versus Genotype/Sex/LifeStage (if they exist in the data\cr
- Mosaic plot of the Sex versus Genotype (if they exist in the data)\cr
}


\seealso{
\code{\link{OpenStatsAnalysis}}, \code{\link{plot.OpenStatsRR}}, \code{\link{plot.OpenStatsMM}}
}

\examples{
####################################################################
# Data preparation
####################################################################
#################
# Categorical data - Creating OpenStatsList object
#################
fileCat <- system.file("extdata", "test_categorical.csv", package = "OpenStats")
test_Cat <- OpenStatsList(
  dataset = read.csv(fileCat, na.strings = "-"),
  testGenotype = "Aff3/Aff3",
  refGenotype = "+/+",
  dataset.colname.genotype = "Genotype",
  dataset.colname.batch = "Assay.Date",
  dataset.colname.lifestage = NULL,
  dataset.colname.weight = "Weight",
  dataset.colname.sex = "Sex"
)
#################
# Fisher's exact test framework
#################
FE_result <- OpenStatsAnalysis(
  OpenStatsList = test_Cat,
  method = "FE",
  FE_formula = Thoracic.Processes ~ Genotype + Sex
)
plot(FE_result, col = 1:2)
}

