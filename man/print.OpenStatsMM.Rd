\name{print.OpenStatsMM}
\alias{print.OpenStatsMM}
\title{
  Summary for an OpenStatsMM object
}
\description{
This function prints summary table for an OpenStatsMM object
}
\usage{
 \method{print}{OpenStatsMM}(x, format = "rst", ...)
}
\arguments{
  \item{x}{
  an instance of OpenStatsMM result from OpenStatsAnalysis(method = 'MM') function
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
	\code{\link{OpenStatsAnalysis}}, 	\code{\link{print.OpenStatsFE}}, 	\code{\link{print.OpenStatsRR}}
}
\examples{
example(OpenStatsAnalysis)
}

