\name{plot.OpenStatsList}
\alias{plot.OpenStatsList}

\title{
	plot for a `OpenStatsList` object
}
\description{
	This function visualises a `OpenStatsList` object
}
\usage{
	\method{plot}{OpenStatsList}(x, vars = NULL, ...)
}

\arguments{
	\item{x}{
		OpenStatsList object 
	}
	\item{vars}{
		Variable(s) of interest. The default is 'Batch', 'Genotype', 'Sex and 'LifeStage' if exists in the data
	}
	\item{\dots}{
		Optional parameters that can be passed to `Hmisc::plot.describe()`
	}
}

\value{Not applicable}

\author{
	Hamed Haseli Mashhadi <hamedhm@ebi.ac.uk>
}

\seealso{
	\code{\link{OpenStatsList}}
}

\examples{
	example(OpenStatsList)
}
