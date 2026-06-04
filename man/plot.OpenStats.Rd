\name{plot.OpenStatsList}
\alias{plot.OpenStatsList}
\title{plot for an `OpenStatsList` object}
\description{This function visualises an `OpenStatsList` object}

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

\details{
The plot function produces two sets of plots for:\cr

- categorical data: scatter plot of proportions\cr
- continuous data: histogram 
}

\value{List of two plot objects, Continuous and Categorical}


\seealso{
	\code{\link{OpenStatsList}}
}

\examples{
example(OpenStatsList)
}
