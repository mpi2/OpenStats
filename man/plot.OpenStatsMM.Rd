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


\value{
Not applicable
}
\details{
The plot function creates some visualisations for the linear mixed model framework: \cr

- Residual versus fitted values \cr 
- Residual density plot and the normality test p-value\cr
- Residual Q-Q plot\cr
- The density plot of the response variable and the normality test p-value\cr
}


\seealso{
\code{\link{OpenStatsAnalysis}}, \code{\link{plot.OpenStatsFE}}, \code{\link{plot.OpenStatsRR}}
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
#  Optimised Linear Mixed model (MM) framework
#################
MM1_result <- OpenStatsAnalysis(
  OpenStatsList = test_Cont,
  method = "MM",
  MM_fixed = data_point ~ Genotype + Weight
)
print(MM1_result, col = 2, main = "Optimised model")
}

