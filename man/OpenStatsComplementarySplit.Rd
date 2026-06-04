\name{OpenStatsComplementarySplit}
\alias{OpenStatsComplementarySplit}
\title{Method "OpenStatsComplementarySplit"}

\description{
This function splits the input data according to the defined values in the `variables` parameter and runs separate analyses on the split datasets. For example, the default split, c("Sex", "LifeStage"), creates independent input data for Males (only), Females, Early, Late, Male.Early, Males.Late, Females.Early, Females.Late and analyses these datasets separately.
}
\usage{
OpenStatsComplementarySplit(
	object    = NULL,
	variables = c("Sex", "LifeStage"),
	debug     = FALSE
)
}
\arguments{
  \item{object}{
  Mandatory argument. An instance of the 'OpenStatsAnalysis' object under the MM (linear mixed model) framework.
}
  \item{variables}{
  Vector of names. A vector of variable names that will be fed into the split engine. The default is `c("Sex", "LifeStage")` that should report the results for the following categories: Males, Females, Early, Late, Male.Early, Males.Late, Females.Early, Females.Late.
}
  \item{debug}{
  Logical flag. Set to TRUE to see the analysis log. Default FALSE.
}
}

\value{
  List of splits and the analysis outputs. The splits contain the name of the partitioning levels (for example Female or Female.Age_15_weeks), and an `OpenStatsAnalysis` object including the input data, outputs etc. See the examples for a general view of the output object.
}



\seealso{
\code{\link{OpenStatsAnalysis}}, \code{\link{plot.OpenStatsComplementarySplit}}, \code{\link{summary.OpenStatsComplementarySplit}}
}
\examples{
####################################################################
# Data preparation
# - Continuous data - Creating OpenStatsList object
####################################################################
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
####################################################################
# Analysis
# -  Optimised Linear Mixed Model (MM) framework
####################################################################
MM_result <- OpenStatsAnalysis(
  OpenStatsList = test_Cont,
  method = "MM",
  MM_fixed = data_point ~ Genotype + Weight
)
####################################################################
# Split on Sex
# ** This split is already available from the normal running of OpenStatsAnalysis
####################################################################
output <- OpenStatsComplementarySplit(object = MM_result, variables = "Sex")
# Structure of the output object
lapply(output, names)
# Summaries
summary(output, format = "pandoc") # See knitr:kable function for more formats
# Plots
plot(output, ask = TRUE)
}
\keyword{ ~OpenStats }
\keyword{ ~SplitEffect }
