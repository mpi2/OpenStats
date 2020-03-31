\name{summary.OpenStatsList}
\alias{summary.OpenStatsList}

\title{
	Summary for a OpenStatsList object
}
\description{
	This function provides a detailed summary of a OpenStatsList object
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

\author{
	Hamed Haseli Mashhadi <hamedhm@ebi.ac.uk>
}

\seealso{
	\code{\link{OpenStatsList}}
}

\examples{
	example(OpenStatsList)
}
