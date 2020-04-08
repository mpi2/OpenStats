\name{OpenStatsListBuilder}
\alias{OpenStatsListBuilder}
\title{
	Method "OpenStatsListBuilder"
}
\description{
	Specifying the age in days, this function creates a 'OpenStatsList' object from a 'PhenList' object. See 'PhenList' function in 'PhenStat' package.
}
\usage{
	OpenStatsListBuilder(
		PhenListobject        ,
		DOE = NULL            ,
		DOB = NULL            ,
		d.threshold = 16 * 7  ,
		debug = TRUE
	)
}
\arguments{
  \item{PhenListobject}{
		'mandatory argument'. Instance of the `PhenList` class
}
  \item{DOE}{
  Name of the data column for the 'Batch' in the `PhenList` object. If left NULL then the input `PhenList` object will be returned. Default NULL
}
  \item{DOB}{
  Name of the data column for the 'date_of_birth' in the `PhenList` object. If left NULL then the input `PhenList` object will be returned. Default NULL
}
  \item{d.threshold}{
  The threshold in age (DOE-DOB)  to specify LifeStage early/late levels. The function uses as.Date(DOE)-as.Date(DOB) to calculate the age. The default is 16 weeks (16*7 days)
}
  \item{debug}{
  Logical flag. Set to TRUE to see debug messages. Default TRUE
}
}

\value{
	Provided DOE and DOB are not NULL, a `OpenStatsList` object that is quite similar to `PhenList` object with an extra column called `LifeStage` with two levels `Early` and `Late`. Otherwise, the output is similar to the input `PhenList` object.
}

\author{
Hamed Haseli Mashhadi <hamedhm@ebi.ac.uk>
}

\seealso{
	\code{\link{OpenStatsAnalysis}}
}
\examples{
	\dontrun{
		file = system.file("extdata", "test_continuous.csv", package = "OpenStats")
		####################################################################
		# PhenListObject
		####################################################################
		PhenListObject = PhenList(
			dataset = read.csv(file),
			testGenotype = 'experimental',
			refGenotype = 'control',
			dataset.colname.batch = 'date_of_experiment',
			dataset.colname.genotype = 'biological_sample_group',
			dataset.colname.sex = 'sex',
			dataset.colname.weight = 'weight'
		)
		####################################################################
		# OpenStatsList object
		####################################################################
		OpenStatsListBuilder   = OpenStatsListBuilder(
			PhenListobject = PhenListObject,
			DOE = 'Batch',
			DOB = 'Birth.Date',
			d.threshold = 99
		)
		plot(OpenStatsListBuilder)
		class(OpenStatsListBuilder)
		rm(OpenStatsListBuilder)
	}
}

\keyword{ OpenStatsListBuilder }

