\name{print.OpenStatsRR}
\alias{print.OpenStatsRR}
\title{
  Summary for an OpenStatsRR object
}
\description{
This function prints summary table for an OpenStatsRR object
}
\usage{
 \method{print}{OpenStatsRR}(x, format = "rst", ...)
}
\arguments{
  \item{x}{
  an instance of OpenStatsRR result from OpenStatsAnalysis(method = 'RR') function
}
  \item{format}{
  See format argument from the knitr::kable function
}
  \item{\dots}{
  Other parameters that can be passed to knitr::kable function
}
}
\value{
The output consists of a pair of values separated by comma, e.g. 1,1, for low and high classes respectively. The following statistics are reported in the summary:

- Applied model \cr
- Checked/optimised model\cr
- Treatment group\cr
- Control group\cr
- If possible, whether sexual dimorphism is detected from the analysis\cr
- Genotype effect p-value\cr
- Genotype effect p-value for females\cr
- Genotype effect p-value for males\cr
- If LifeStage existed in the data, LifeStage p-value\cr
- Genotype effect for early adults\cr
- Genotype effect for late adults\cr
- If Sex existed in the data, Sex p-value\cr
- If bodyweight existed in the data, bodyweight p-value\cr
}


\seealso{
	\code{\link{OpenStatsAnalysis}},, 	\code{\link{print.OpenStatsFE}}, 	\code{\link{print.OpenStatsMM}}
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
print(RR_result)
}

