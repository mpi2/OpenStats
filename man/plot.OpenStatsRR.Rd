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

\value{
Not applicable
}

\details{
The plot function creates some visualisations for the reference range plus framework \cr

- Mosaic plot of the discretised response versus Genotype/Sex/LifeStage (if they exist in the data\cr
- Mosaic plot of the Sex versus Genotype (if they exist in the data)\cr
}

\seealso{
\code{\link{OpenStatsAnalysis}}, \code{\link{plot.OpenStatsFE}}, \code{\link{plot.OpenStatsMM}}
}

\examples{
####################################################################
# Data preparation
####################################################################
#################
# Continuous data - Creating OpenStatsList object
#################
fileCon <- system.file("extdata", "test_continuous.csv", package = "OpenStats")
test_Cont <- OpenStatsList(
  dataset = read.csv(fileCon),
  testGenotype = "experimental",
  refGenotype = "control",
  dataset.colname.genotype = "biological_sample_group",
  dataset.colname.batch = "date_of_experiment",
  dataset.colname.lifestage = NULL,
  dataset.colname.weight = "weight",
  dataset.colname.sex = "sex"
)

#################
# Reference range framework
#################
RR_result <- OpenStatsAnalysis(
  OpenStatsList = test_Cont,
  method = "RR",
  RR_formula = data_point ~ Genotype + Sex
)
plot(RR_result, col = 3:4)
}

