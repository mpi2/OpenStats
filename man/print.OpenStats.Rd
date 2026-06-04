\name{print.OpenStatsList}
\alias{print.OpenStatsList}

\title{
Print summary table for an OpenStatsList object
}
\description{
	This function prints a summary table for an OpenStatsList object
}
\usage{
	\method{print}{OpenStatsList}(x, vars = NULL, ...)
}


\arguments{
	\item{x}{
		OpenStatsList object 
	}
	\item{vars}{
		Variable(s) of interest
	}
	\item{\dots}{
		Optional parameters that can be passed to 'summarytools::dfSummary ()'
	}
}
\value{
Tabe of summary statistics
}


\seealso{
	\code{\link{OpenStatsList}}, \code{\link{summary.OpenStatsList}}, \code{\link{OpenStatsAnalysis}}
}

\examples{
example(OpenStatsList)
}
