\name{OpenStatsAnalysis}
\alias{OpenStatsAnalysis}
\title{Method "OpenStatsAnalysis"}

\description{
  The driver function in the OpenStats package for running statisical analysis on phenotypic data. \cr

- It performs several checks on the data and the input model before performing the analysis. This function supports three main analysis frameworks precisely, Linear Mix model (MM), Fisher's Exact test (FE) and Reference Range plus (RR).\cr
- It further monitors the process for failures and errors and applies some runtime patches/fixes.\cr
- The function parameters are designed to be human-friendly by initialising the inputs by the model that will be applied to the data.

}
\usage{
OpenStatsAnalysis(
	OpenStatsListObject    = NULL,
	method     = NULL,
	MM_fixed   = TypicalModel(
		depVariable = "data_point",
		withWeight  = MM_BodyWeightIncluded,
		Sex = TRUE,
		LifeStage = TRUE,
		data = OpenStatsListObject@datasetPL,
		others = NULL,
		debug = debug
	),
	MM_random = rndProce("TYPICAL"),
	MM_BodyWeightIncluded = TRUE,
	MM_lower  = ~ Genotype + 1,
	MM_weight = if (
		TermInModelAndnLevels(
				 model = MM_fixed,
				 data = OpenStatsListObject@datasetPL
				 )
			){
		varIdent(form = ~ 1 |  LifeStage)
	}else{
		varIdent(form = ~ 1 |  Genotype)
		},
	MM_direction = "both",
	MM_checks    = c(TRUE, TRUE, TRUE, TRUE),
	MM_optimise  = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
	FE_formula   = category ~   Genotype + Sex + LifeStage,
	RR_formula   = data_point ~ Genotype + Sex + LifeStage,
	RRrefLevel = 'control',
	RR_prop      = 0.95,
	FERR_rep     = 1500,
	FERR_FullComparisions = c(TRUE, FALSE),
	MMFERR_conf.level = 0.95,
	debug        = TRUE,
	...
)
}
\arguments{
\item{OpenStatsListObject}{
mandatory argument. An instance of the 'OpenStatsList' or 'PhenList' (former from the Bioconductor `PhenStat` package).
}
\item{method}{
Must be specified by the user. A character string ("MM", "FE" or "RR") defining the method to use for model building.
}
\item{MM_fixed}{
Only applies to the "MM" framework. A formula that specifies the fixed effect in the linear mix model. The default is data_point~Genotype+Sex+LifeStage +/- BodyWeight. Note that the algorithm checks and removes the formula terms that do not exist in data.
}
\item{MM_random}{
Only applies to the "MM" framework. The random effect in the linear mixed model. See the lme() function for more details about the formula. The default is ~1|Batch.
}
\item{MM_BodyWeightIncluded}{
Only applies to the `default` MM_fixed in the "MM" framework. If TRUE then the default model includes the body weight and the model would be: data_point~Genotype+Sex+LifeStage + Weight.
}
\item{MM_lower}{
Only applies to the "MM" framework. A right-sided formula, for example ~Genotype+1 or ~Sex+Genotype+1 or ~Sex+Genotype+Sex:Genotype. The lowest model that must not be included in the model optimisation. In other words, the terms in this model won't be removed during the optimisation process. The default is ~ Genotype + 1 that is the genotype effect and the intercept will be kept in the model during the optimisation process.
}
\item{MM_weight}{
Only applies to the "MM" framework. From weight in the lme() manual: \cr
"an optional varFunc object or one-sided formula describing the within-group heteroscedasticity structure. If given as a formula, it is used as the argument to varFixed, corresponding to fixed variance weights. See the documentation on varClasses for a description of the available varFunc classes. Defaults to NULL, corresponding to homoscedastic within-group errors". \cr\cr
The default is varIdent(form = ~ 1 |  LifeStage) if the LifeStage included in the input data. Otherwise,  varIdent(form = ~ 1 | Genotype).
}
\item{MM_direction}{
Only applies to the "MM" framework. Select from  "both" (for stepwise optimisation), "backward" (for backward elimination) or "forward" (for forward selection) for the optimisation algorithm. The default is "both".
}
\item{MM_checks}{
Only applies to the "MM" framework. A vector of four 1/0 or TRUE/FALSE values such as c(TRUE, TRUE, TRUE, TRUE)[default]. Performing pre checks on the input model for some known scenarios. The first element of the vector activates checks on the model terms (See MM_fixed) to be existed in data. The second term removes any single level -factor- from the model (in MM_fixed). The third term removes the single value (such as a column of constants/no variation) from the -continuous- terms in the model (in MM_fixed). The Fourth element checks the interaction term to make sure all interactions have some data attached. Caution is needed for this check as it may take longer than usual if the formula in MM_fixed contains many factors. The default is c(TRUE, TRUE, TRUE, TRUE) that is all checks perform.\cr

* Note that the function always removes duplicated columns in the dataset prior to applying the lme/gls.\cr
* Regardless of the `check` settings, the function always checks for the existence of the `MM_random` terms (given it is not set to NULL) in the input data
}
\item{MM_optimise}{
Only applies to the "MM" framework. A vector of six binary values such as c(1,1,1,1,1,1) or c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) (default). The first element of the vector activates the fixed effect optimisation. The algorithm uses AICc to optimise the fixed effects (Check `AICcmodavg` package for more details about AICc). The second and third elements of the vector activate optimisation on 'weight' and 'random effects' respectively. The optimisation of weight and random effects refers to comparing the AICc between a model with and without those effects. The fourth element activates the Split model effects (for example, separate male and female effects) (see `SplitModels` in the output object). The fifth effect activates the effect size estimation (see `Effect sizes` in the output object). The sixth element activates the normality tests on the residuals (see `ResidualNormalityTests` in the output object).
}
\item{FE_formula}{
Only applies to the "FE" framework. The model for analysing the categorical data. The default is: category ~ Genotype + Sex + LifeStage. \cr
Note that similar to MM_fixed, the terms that do not exist in data or the interaction terms will be dismissed from the model.
}
\item{RR_formula}{
Only applies to the "RR" framework. The model for analysing the RR+ compatible data. \cr
** Important. The first term on the right hand side of the formula specifies the variable for discritising the response.\cr
The default is: data_point ~ Genotype + Sex + LifeStage. \cr
Note that similar to MM_fixed,  the terms that do not exist in data or the interaction terms will be dismissed from the model.
}
\item{RRrefLevel}{
Only applies to the "RR" framework. A single term for the `reference level` in the `reference variable` (the first term on the right hand side of the `RR_formula`) used for discritising the response. If left blank then the a level with more observations will be considered as the reference level. The default is `control`.
}  \item{RR_prop}{
Only applies to the "RR" framework. A single value between (0.5,1) not including the boundaries. The threshold for the variation ranges in the RR framework. The default value is 0.95.
}
\item{FERR_rep}{
Only applies to the "RR" or "FE" frameworks. The number of iteration for the Monte Carlo Fisher's Exact test. See "B" parameter in `fisher.test()` function. Set to 0 for non-bayesian results (not recommended). The default is 1500.
}
\item{FERR_FullComparisions}{
Only applies to the "RR" or "FE" frameworks. A vector of two logical flags, default c(TRUE,FALSE). Setting the first value to TRUE, then all combinations of the effects (all levels of factors in the input model - for example Male_LifeStage, Male_Genotype, Male_Mutant, Male_control, Female_control, Female_Mutant, Female_LifeStage and so on) will be tested. Otherwise only main effects (no sub levels - for example Sex_LifeStage [not for instance Male_LifeStage]) will be tested. Setting the second element of the vector to TRUE (default FALSE) will force Fisher's Exact test to do all comparisions between different levels of the RESPONSE variable. For example, if the respose has three levels such as 1.positive, 2.negative and 3.neutral then the comparisions will be between  1\&2, 1\&3, 2\&3 and 1\&2\&3 (obviously this is the full table).
}
\item{MMFERR_conf.level}{
Applies to all frameworks (MM, FE, RR). Single numeric value for the interval confidence level. Default is 0.95
}
\item{debug}{
A logical flag. Set to TRUE to see more details about the progress of the function. Default TRUE
}
\item{\dots}{
Other parameters that can be passed to:\cr
~> If the model is set to Linear Mixed Model (MM) then the parameters that can be passed to `lme` function. See ?lme() manual page\cr
~> If the model is either Fisher's Exact test (FE) or Reference Range + (RR) then the parameters that can be passed to the `fisher.test` function. See ?fisher.test() manual page\cr
}
}

\value{
  1. Successful execution of the function will return a list of three elements:
  \item{input }{This contains the list of inputs}
  \item{output }{A list of outputs}
  \item{extra }{A placeholder for extra information if exists}
  2. If the function fails:
  \item{messages}{A placeholder for the errors/warnings in the case of failure}
}

\details{
OpenStatsReport function can be used to extract the key elements of the analysis from the OpenStatsMM/FE/RR objects. The output from OpenStatsReport has schemed that makes it easy to be populated to the downstream processes such as storing and accessing results from a database.
}
\seealso{
  \code{\link{OpenStatsList}}\code{\link{OpenStatsComplementarySplit}}, \code{\link{plot.OpenStatsMM}},\code{\link{plot.OpenStatsFE}},\code{\link{plot.OpenStatsRR}}, \code{\link{summary.OpenStatsMM}}, \code{\link{summary.OpenStatsFE}},\code{\link{summary.OpenStatsRR}}, \code{\link{print.OpenStatsMM}}, \code{\link{print.OpenStatsFE}}, \code{\link{print.OpenStatsRR}}
}
\examples{
####################################################################
# 1 Data preparation
####################################################################
#################
# 1.1 Continuous data - Creating OpenStatsList object
#################
fileCon <- system.file("extdata", "test_continuous.csv", package = "OpenStats")
test_Cont <- OpenStatsList(
  dataset = read.csv(fileCon),
  testGenotype = "experimental",
  refGenotype = "control",
  dataset.colname.genotype = "biological_sample_group",
  dataset.colname.batch = "date_of_experiment",
  dataset.colname.lifestage = NULL,
  dataset.colname.weight = "weight",
  dataset.colname.sex = "sex"
)
#################
# 1.2 Categorical data - Creating OpenStatsList object
#################
fileCat <- system.file("extdata", "test_categorical.csv", package = "OpenStats")
test_Cat <- OpenStatsList(
  dataset = read.csv(fileCat, na.strings = "-"),
  testGenotype = "Aff3/Aff3",
  refGenotype = "+/+",
  dataset.colname.genotype = "Genotype",
  dataset.colname.batch = "Assay.Date",
  dataset.colname.lifestage = NULL,
  dataset.colname.weight = "Weight",
  dataset.colname.sex = "Sex"
)
####################################################################
# 2 Testing frameworks
####################################################################

#################
# 2.1 Optimised Linear Mixed model (MM) framework
#################
MM1_result <- OpenStatsAnalysis(
  OpenStatsList = test_Cont,
  method = "MM",
  MM_fixed = data_point ~ Genotype + Weight
)
VO_MM1 <- OpenStatsReport(MM1_result)
plot(MM1_result, col = 2, main = "Optimised model")
summary(MM1_result)


#################
# 2.2 Linear Mixed model (MM) with NO optimisation
# for the fixed effects but random/weight effects
#################
MM2_result <- OpenStatsAnalysis(
  OpenStatsList = test_Cont,
  method = "MM",
  MM_fixed = data_point ~ Genotype + Weight,
  MM_lower = ~ Genotype + Weight + 1
  # Or simply MM_optimise = c(0, 1, 1, 1, 1, 1)
)
VO_MM2 <- OpenStatsReport(MM2_result)
plot(MM2_result, col = 8, main = "No optimisation on the fixed effects")
summary(MM2_result)

#################
# 2.3 Linear Mixed model (MM) with NO optimisation on the model
#################
MM3_result <- OpenStatsAnalysis(
  OpenStatsList = test_Cont,
  method = "MM",
  MM_fixed = data_point ~ Genotype + Weight,
  MM_optimise = c(0, 0, 0, 1, 1, 1)
)
VO_MM3 <- OpenStatsReport(MM3_result)
plot(MM3_result, col = 3, main = "Not optimised model")
summary(MM3_result)


#################
# 2.4 Reference range framework
#################
RR_result <- OpenStatsAnalysis(
  OpenStatsList = test_Cont,
  method = "RR",
  RR_formula = data_point ~ Genotype + Sex
)
VO_RR <- OpenStatsReport(RR_result)
plot(RR_result, col = 3:4)
summary(RR_result)


#################
# 2.5 Fisher's exact test framework
#################
FE_result <- OpenStatsAnalysis(
  OpenStatsList = test_Cat,
  method = "FE",
  FE_formula = Thoracic.Processes ~ Genotype + Sex
)
VO_FE <- OpenStatsReport(FE_result)
plot(FE_result, col = 1:2)
summary(FE_result)
}
\keyword{ ~OpenStats }
