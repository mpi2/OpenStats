# testing OpenStats
# Author: Hamed Haseli Mashhadi
# Maintainer: Marina Kan <marinak@ebi.ac.uk>
##########################################################################
################### General Testing ######################################

##########################################################################
# In all tests below we compare OpenStats output with PhenStat
# PhenStat: https://www.bioconductor.org/packages/release/bioc/html/PhenStat.html
##########################################################################


##########################################################################
# Prepare the environment
##########################################################################
library(OpenStats)
# If package PhenStat does exist, then install the package from Bioconductor
if (!"PhenStat" %in% installed.packages()[, "Package"]) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install("PhenStat")
} else {
  library(PhenStat)
}

as.numericRound <- function(x, digit = 8) {
  r <- round(as.numeric(x), digits = digit)
  return(r)
}

##########################################################################
# 1. Basic checks
##########################################################################
test_OpenStatsList <- function() {
  message("Testing the basic checks ...")

  # Load the raw data
  fileCon <- system.file("extdata", "test_continuous.csv", package = "OpenStats")
  data_continuous <- read.csv(fileCon, as.is = TRUE)

  # The data object from OpenStats
  OpenStatsListObject <-
    OpenStatsList(
      dataset = data_continuous,
      testGenotype = "experimental",
      refGenotype = "control",
      dataset.colname.genotype = "biological_sample_group",
      dataset.colname.batch = "date_of_experiment",
      dataset.colname.sex = "sex",
      dataset.colname.weight = "weight",
      debug = FALSE
    )

  # The data object from PhenStat
  PhenStatListObject <-
    PhenList(
      dataset = data_continuous,
      testGenotype = "experimental",
      refGenotype = "control",
      dataset.colname.genotype = "biological_sample_group",
      dataset.colname.batch = "date_of_experiment",
      dataset.colname.sex = "sex",
      dataset.colname.weight = "weight",
      outputMessages = FALSE
    )

  compare <- identical(
    OpenStatsListObject@datasetPL,
    PhenStatListObject@datasetPL
  )
  if (!compare) {
    return("test failed")
  } else {
    return("test passed")
  }
}


############################################################################
# 2. Test Reference Range plus framework
############################################################################
test_RR_Framework <- function() {
  message("Testing the Reference Range plus framework ...")

  # Load the raw data
  fileCon <- system.file("extdata", "test_continuous.csv", package = "OpenStats")
  data_continuous <- read.csv(fileCon, as.is = TRUE)

  # The data object from OpenStats
  OpenStatsListObject <-
    OpenStatsList(
      dataset = data_continuous,
      testGenotype = "experimental",
      refGenotype = "control",
      dataset.colname.genotype = "biological_sample_group",
      dataset.colname.batch = "date_of_experiment",
      dataset.colname.sex = "sex",
      dataset.colname.weight = "weight",
      debug = FALSE
    )

  # The data object from PhenStat
  PhenStatListObject <-
    PhenList(
      dataset = data_continuous,
      testGenotype = "experimental",
      refGenotype = "control",
      dataset.colname.genotype = "biological_sample_group",
      dataset.colname.batch = "date_of_experiment",
      dataset.colname.sex = "sex",
      dataset.colname.weight = "weight",
      outputMessages = FALSE
    )

  # The analysis results from OpenStats
  testOpenStats <- OpenStatsAnalysis(
    OpenStatsListObject = OpenStatsListObject,
    method = "RR",
    debug = FALSE
  )

  # The analysis results from PhenStat
  testPhenStat <- testDataset(
    phenList = PhenStatListObject,
    depVariable = "data_point",
    equation = "withoutWeight",
    method = "RR",
    outputMessages = FALSE
  )

  # Vectorising the outputs
  VectoriseOutputPhenStat <- as.list(vectorOutput(testPhenStat))
  VectoriseOutputOpenStats <- OpenStatsReport(testOpenStats)

  # Comparing the outputs
  compare <- identical(
    VectoriseOutputPhenStat$`Genotype p-Val`,
    paste(
      as.numericRound(VectoriseOutputOpenStats$`Genotype p-value`$Low$p.value),
      VectoriseOutputOpenStats$`Genotype p-value`$High$p.value,
      sep = ","
    )
  )

  if (!compare) {
    return("test failed")
  } else {
    return("test passed")
  }
}


############################################################################
# 3. Test Fisher's exact test framework
############################################################################
test_FE_Framework <- function() {
  message("Testing Fisher's Exact test framework ...")

  # Load the raw data
  fileCat <- system.file("extdata", "test_categorical.csv", package = "OpenStats")
  data_categorical <-
    read.csv(fileCat, na.strings = "-", as.is = TRUE)

  # The data object from OpenStats
  OpenStatsListObject <-
    OpenStatsList(
      dataset = data_categorical,
      testGenotype = "Aff3/Aff3",
      refGenotype = "+/+",
      dataset.colname.genotype = "Genotype",
      dataset.colname.batch = "Assay.Date",
      dataset.colname.weight = "Weight",
      dataset.colname.sex = "Sex",
      debug = FALSE
    )

  # The data object from PhenStat
  PhenStatListObject <-
    PhenList(
      dataset = data_categorical,
      testGenotype = "Aff3/Aff3",
      refGenotype = "+/+",
      dataset.colname.genotype = "Genotype",
      dataset.colname.batch = "Assay.Date",
      dataset.colname.weight = "Weight",
      dataset.colname.sex = "Sex",
      outputMessages = FALSE
    )

  # The analysis results from OpenStats
  testOpenStats <- OpenStatsAnalysis(
    OpenStatsListObject = OpenStatsListObject,
    FE_formula = Thoracic.Processes ~ Genotype + Sex,
    method = "FE",
    FERR_rep = 0,
    debug = FALSE
  )

  # The analysis results from PhenStat
  testPhenStat <- testDataset(
    phenList = PhenStatListObject,
    depVariable = "Thoracic.Processes",
    equation = "withoutWeight",
    method = "RR",
    outputMessages = FALSE
  )

  # Comparing the outputs
  compare <- identical(
    testOpenStats$output$SplitModels$Thoracic.Processes$Genotype$result$`Complete table`$p.value,
    testPhenStat@analysisResults[[1]]@modelOutput$p.value
  )

  if (!compare) {
    return("test failed")
  } else {
    return("test passed")
  }
}



############################################################################
# 4. Test Linear mixed model framework
############################################################################
test_MM_Framework <- function() {
  message("Testing Linear Mixed model framework ...")

  # Load the raw data
  fileCon <- system.file("extdata", "test_continuous.csv", package = "OpenStats")
  data_continuous <- read.csv(fileCon, as.is = TRUE)

  # The data object from OpenStats
  OpenStatsListObject <-
    OpenStatsList(
      dataset = data_continuous,
      testGenotype = "experimental",
      refGenotype = "control",
      dataset.colname.genotype = "biological_sample_group",
      dataset.colname.batch = "date_of_experiment",
      dataset.colname.sex = "sex",
      dataset.colname.weight = "weight",
      debug = FALSE
    )

  # The data object from PhenStat
  PhenStatListObject <-
    PhenList(
      dataset = data_continuous,
      testGenotype = "experimental",
      refGenotype = "control",
      dataset.colname.genotype = "biological_sample_group",
      dataset.colname.batch = "date_of_experiment",
      dataset.colname.sex = "sex",
      dataset.colname.weight = "weight",
      outputMessages = FALSE
    )

  # The analysis results from OpenStats
  testOpenStats <- OpenStatsAnalysis(
    OpenStatsListObject = OpenStatsListObject,
    method = "MM",
    data_point ~ Genotype + Sex + Weight,
    MM_BodyWeightIncluded = TRUE,
    MM_optimise = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
    debug = FALSE
  )

  # The analysis results from PhenStat
  testPhenStat <- testDataset(
    phenList = PhenStatListObject,
    depVariable = "data_point",
    equation = "withWeight",
    method = "MM",
    outputMessages = FALSE,
    keepList = c(
      keep_batch = TRUE,
      keep_equalvar = TRUE,
      keep_weight = TRUE,
      keep_sex = TRUE,
      keep_interaction = FALSE
    )
  )

  # Vectorising the outputs
  VectoriseOutputOpenStats <- OpenStatsReport(testOpenStats)
  VectoriseOutputPhenStat <- as.list(vectorOutput(testPhenStat))

  # Comparing the Genotype effect p-values
  VectoriseOutputPhenStat$`Genotype p-Val`
  VectoriseOutputOpenStats$`Genotype contribution`$Overall

  compare <- abs(
    VectoriseOutputOpenStats$`Genotype p-value` - as.numericRound(VectoriseOutputPhenStat$`Genotype p-Val`)
  ) <
    0.005

  # Comparing the outputs
  if (!compare) {
    return("test failed")
  } else {
    cat(
      "Test passed with toleration = 0.005. \nJustification:\n - The increase in the number of iterations in OpenStats (1500) versus PhenStat (50)\n."
    )
    return("test passed")
  }
}

test_OpenStatsList()
test_RR_Framework()
test_FE_Framework()
test_MM_Framework()
