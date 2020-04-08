\name{OpenStatsList}
\alias{OpenStatsList}

\title{
Method "OpenStatsList"
}
\description{
Function to create `OpenStatsList` object from a data frame.

~> The mandatory variables for creating a `standard` OpenStatsList objects are `response variable` and `Genotype`. Having two levels in the `Genotype` is mandatory. The function further checks for the `Sex` with two levels (Male/Female), LifeStage` with two levels (Early/Late), `Batch` (defined as date_of experiment in the IMPC) and `Weight` (defined as mice body weight in the IMPC) and reports any abnormality in the data. \cr

~> For advance applications, the function is capable of creating a `OpenStatsList` object without performing  checks. To do this, set clean.dataset = FALSE.\cr
}
\usage{
OpenStatsList(
	dataset                                                 ,
	testGenotype                = 'experimental'            ,
	refGenotype                 = 'control'                 ,
	hemiGenotype                = NULL                      ,
	clean.dataset               = TRUE                      ,
	dataset.colname.genotype    = 'biological_sample_group' ,
	dataset.colname.sex         = 'sex'                     ,
	dataset.colname.batch       = 'date_of_experiment'      ,
	dataset.colname.lifestage   = 'LifeStage'               ,
	dataset.colname.weight      = 'weight'                  ,
	dataset.values.missingValue = c(' ','')                 ,
	dataset.values.male         = NULL                      ,
	dataset.values.female       = NULL                      ,
	dataset.values.early        = NULL                      ,
	dataset.values.late         = NULL                      ,
	debug                       = TRUE
)
}

\arguments{
  \item{dataset}{
  'mandatory argument'. data frame created from file or from another source
}
  \item{testGenotype}{
  `mandatory argument`. Defines the test genotype to be compared to the reference genotype. Default `experimental`
}
  \item{refGenotype}{
  defines the reference genotype; assigned default value is `control`
}
  \item{hemiGenotype}{
  defines the genotype value for hemizygous that will be changed to test genotype value
}
  \item{clean.dataset}{
  logical flag. `TRUE` activates all checks and modification on the input data. The overview of the checks is, existence of the variables, checking levels, missings and relabeling
}
 \item{dataset.colname.genotype}{
  column name within dataset for the genotype. Default 'biological_sample_group'
}
  \item{dataset.colname.sex}{
  column name within dataset for the sex. Default 'sex'
}
  \item{dataset.colname.batch}{
  column name within dataset for the batch effect. Default 'date_of_experiment'
}
  \item{dataset.colname.lifestage}{
  column name within dataset for the life stage. Default 'LifeStage'
}
  \item{dataset.colname.weight}{
  column name within dataset for the body weight. Default 'weight'
}
  \item{dataset.values.missingValue}{
  value used as missing value in the dataset. Default `(space)`.
}
  \item{dataset.values.male}{
  value used to label "males" in the dataset
}
  \item{dataset.values.female}{
  value used to label "females" in the dataset
}
  \item{dataset.values.early}{
  value used to label "early life stage" in the dataset
}
  \item{dataset.values.late}{
  value used to label "late life stage" in the dataset
}
  \item{debug}{
	A logical flag. Set to TRUE to see more details about the progress of the function. Default TRUE
}
}
\value{
an instance of the OpenStatsList class.
}

\author{
Hamed Haseli Mashhadi <hamedhm@ebi.ac.uk>
}
\seealso{
	\code{\link{OpenStatsAnalysis}}
}
\examples{
  ####################################################################
  df = read.csv(system.file("extdata", "test_continuous.csv", package = "OpenStats"))
  ####################################################################
  # OpenStatsList object
  ####################################################################
  OpenStatsList   = OpenStatsList(
  	dataset   = df,
  	testGenotype = 'experimental',
  	refGenotype = 'control',
  	dataset.colname.batch = 'date_of_experiment',
  	dataset.colname.genotype = 'biological_sample_group',
  	dataset.colname.sex = 'sex',
  	dataset.colname.weight = 'weight'
  )
  plot   (OpenStatsList)
  summary(OpenStatsList, style = 'grid')
  class  (OpenStatsList)
  rm     (OpenStatsList)
  # str(OpenStatsList@datasetPL)
}

\keyword{ OpenStatsList }
\keyword{ OpenStatsAnalysis }
