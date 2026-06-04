\name{summary.OpenStatsList}
\alias{summary.OpenStatsList}

\title{
	Summary for an OpenStatsList object
}
\description{
	This function provides a detailed summary of an OpenStatsList object
}
\usage{
	\method{summary}{OpenStatsList}(object, vars = NULL, ...)
}


\arguments{
	\item{object}{
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
	\code{\link{OpenStatsList}}, \code{\link{summary.OpenStatsList}}
}

\examples{
example(OpenStatsList)
}
