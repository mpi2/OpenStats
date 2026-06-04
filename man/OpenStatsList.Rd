\name{OpenStatsList}
\alias{OpenStatsList}

\title{Method "OpenStatsList"}

\description{
The driver function to create `OpenStatsList` object from a data frame.

- The mandatory variable for creating a `standard` OpenStatsList objects is `Genotype`. Having two levels in the `Genotype` field is mandatory. The function further checks for the optional `Sex` with two levels (Male/Female), LifeStage` with two levels (Early/Late), `Batch` (defined as date_of experiment in the IMPC) and `Weight` (defined as animal body weight in the IMPC) and reports any abnormality in the data. \cr

- For advance applications, the function is capable of creating a `OpenStatsList` object without performing  checks. To do this, set clean.dataset to FALSE.\cr
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
  mandatory argument. data frame created from file or from another source. See notes for more details
}
  \item{testGenotype}{
  mandatory argument. Defines the test genotype to be compared to the reference genotype. Default `experimental`
}
  \item{refGenotype}{
  defines the reference genotype; assigned default value is `control`
}
  \item{hemiGenotype}{
  optional argument. defines the genotype value for hemizygous that will be changed to test genotype value
}
  \item{clean.dataset}{
  logical flag. `TRUE` activates all checks and modification on the input data. The overview of the checks is, existence of the variables, checking levels, missings and relabeling
}
 \item{dataset.colname.genotype}{
  mandatory argument. Column name within dataset for the genotype. Default 'biological_sample_group'
}
  \item{dataset.colname.sex}{
  optional argument. column name within dataset for the sex. Default 'sex'
}
  \item{dataset.colname.batch}{
  optional argument. column name within dataset for the batch effect. Default 'date_of_experiment'
}
  \item{dataset.colname.lifestage}{
  optional argument. column name within dataset for the life stage. Default 'LifeStage'
}
  \item{dataset.colname.weight}{
  optional argument. column name within dataset for the body weight. Default 'weight'
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

\note{
OpenStats allows a `data.frame` for the input data. This data.frame can be formed from csv, tsv, txt etc. files and is organised with rows and columns for samples and features respectively.  This allows a wide range of integration with other Bioconductor/CRAN packages, for instance, the output of Bioconductor `SummarizedExperiment` package can be transformed and fed into OpenStats (note that SummarizedExperiment allows sample in columns and feature in rows that requires at least a transpose operation). Additionally, Bioconductor `PhenStat` function `PhenList` produces very similar results to `OpenStatsList` that allows direct processing of the `PhenList` object by downstream OpenStats operational functions. 
}

\value{
an instance of the OpenStatsList class. The S4 object contains:\cr

1. raw data: `OpenStatsListObject@datasetUNF`\cr
2. polished `data: OpenStatsListObject@datasetPL`\cr
3. the inputarguments to the `OpenStatsList` function\cr
}

\seealso{
	\code{\link{OpenStatsAnalysis}},\code{\link{plot.OpenStatsList}},\code{\link{summary.OpenStatsList}}, \code{\link{summary.OpenStatsList}},
}
\examples{
####################################################################
df <- read.csv(system.file("extdata", "test_continuous.csv", package = "OpenStats"))
####################################################################
# OpenStatsList object
####################################################################
OpenStatsList <- OpenStatsList(
  dataset = df,
  testGenotype = "experimental",
  refGenotype = "control",
  dataset.colname.batch = "date_of_experiment",
  dataset.colname.genotype = "biological_sample_group",
  dataset.colname.sex = "sex",
  dataset.colname.weight = "weight"
)
p <- plot(OpenStatsList,
  vars = c(
    "Genotype",
    "Sex",
    "data_point",
    "age_in_days"
  )
)
p$Continuous
p$Categorical
summary(OpenStatsList, style = "grid")
class(OpenStatsList)
rm(OpenStatsList)
}

\keyword{ OpenStatsList }
\keyword{ OpenStatsAnalysis }
