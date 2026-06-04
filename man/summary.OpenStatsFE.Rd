\name{summary.OpenStatsFE}
\alias{summary.OpenStatsFE}
\title{
  Summary for an OpenStatsFE object
}
\description{
This function provides summary for an OpenStatsFE object
}
\usage{
 \method{summary}{OpenStatsFE}(object, format = "rst", ...)
}
\arguments{
  \item{object}{
  an instance of OpenStatsFE result from OpenStatsAnalysis(method = 'FE') function
}
  \item{format}{
  See format argument from the knitr::kable function
}
  \item{\dots}{
  Other parameters that can be passed to knitr::kable function
}
}

\value{
The output consists of the following statistics:

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
	\code{\link{OpenStatsAnalysis}}, 	\code{\link{summary.OpenStatsMM}}, 	\code{\link{summary.OpenStatsRR}}
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
summary(FE_result)
}

