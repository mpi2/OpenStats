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

\author{
	Hamed Haseli Mashhadi <hamedhm@ebi.ac.uk>
}

\seealso{
	\code{\link{OpenStatsAnalysis}}, 	\code{\link{summary.OpenStatsMM}}, 	\code{\link{summary.OpenStatsRR}}
}
\examples{
example(OpenStatsAnalysis)
}

