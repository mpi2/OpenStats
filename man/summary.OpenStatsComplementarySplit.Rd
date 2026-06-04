\name{summary.OpenStatsComplementarySplit}
\alias{summary.OpenStatsComplementarySplit}
\title{
  Summary for an OpenStatsComplementarySplit object
}
\description{
This function provides summary for an `OpenStatsComplementarySplit` object
}
\usage{
 \method{summary}{OpenStatsComplementarySplit}(object, format = "rst", ...)
}
\arguments{
  \item{object}{
  an instance of `OpenStatsComplementarySplit` result 
}
  \item{format}{
  See format argument from the knitr::kable function
}
  \item{\dots}{
  Other parameters that can be passed to knitr::kable function
}
}
\value{
The output consists of the following statistics for levels of partitioning variables (see `variables` in the `OpenStatsComplementarySplit` function manual):

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
	\code{\link{OpenStatsComplementarySplit}}, 	\code{\link{plot.OpenStatsComplementarySplit}}, \code{\link{print.OpenStatsComplementarySplit}}, \code{\link{OpenStatsAnalysis}}
}
\examples{
example(OpenStatsComplementarySplit)
}

