OpenStatsReportRR <- function(object) {
  if (!is.null(object$messages)) {
    return(NULL)
  }
  #####################################################################
  Labels <- OpenStatsListLevels(object = object)
  Fmodel <- object$extra$Cleanedformula
  frm <- formula(Fmodel)
  depVariable <- all_vars0(frm)[1]
  refVariable <- all_vars0(frm)[2]
  # 	equation     = NULL
  formula <- printformula(frm)
  framework <- paste0(
    "Reference Range Plus Test framework; quantile = ",
    object$input$prop,
    " (Tails probability = ",
    1 - MakeRRQuantileFromTheValue(object$input$prop, messages = FALSE),
    ")"
  )
  # 	fittingMethod  = NULL
  #####################################################################
  Vsplit <- object$output$SplitModels
  low <- pastedot("Low", depVariable, refVariable)
  high <- pastedot("High", depVariable, refVariable)
  VsplitLow <- Vsplit[[low]]
  VsplitHig <- Vsplit[[high]]
  #### Sex
  VsplitLowFemale <- Vsplit[[pastedot(low, "Female")]]
  VsplitHigFemale <- Vsplit[[pastedot(high, "Female")]]
  VsplitLowMale <- Vsplit[[pastedot(low, "Male")]]
  VsplitHigMale <- Vsplit[[pastedot(high, "Male")]]
  ### Lifestage
  VsplitLowEarly <- Vsplit[[pastedot(low, "Early")]]
  VsplitHigEarly <- Vsplit[[pastedot(high, "Early")]]
  VsplitLowLate <- Vsplit[[pastedot(low, "Late")]]
  VsplitHigLate <- Vsplit[[pastedot(high, "Late")]]
  ### LifeStage x Sex
  # Female
  VsplitLowEarlyFemale <- Vsplit[[pastedot(low, "Female", "Early")]]
  VsplitHigEarlyFemale <- Vsplit[[pastedot(high, "Female", "Early")]]
  VsplitLowLateFemale <- Vsplit[[pastedot(low, "Female", "Late")]]
  VsplitHigLateFemale <- Vsplit[[pastedot(high, "Female", "Late")]]
  # Male
  VsplitLowEarlyMale <- Vsplit[[pastedot(low, "Male", "Early")]]
  VsplitHigEarlyMale <- Vsplit[[pastedot(high, "Male", "Early")]]
  VsplitLowLateMale <- Vsplit[[pastedot(low, "Male", "Late")]]
  VsplitHigLateMale <- Vsplit[[pastedot(high, "Male", "Late")]]
  ###
  GenotypeDiscLabel <- "Data is discritised by `Genotype` levels"
  SexDiscLabel <- "Data is first split by `Sex` levels then discritised by `Genotype` levels"
  LifeStageDiscLabel <- "Data is first split by `LifeStage` levels then discritised by `Genotype` levels"
  LifeStageSexDiscLabel <- "Data is first split by `LifeStage-Sex` levels then discritised by `Genotype` levels"
  #####################################################################
  x <- object$input$data
  columnOfInterest <- x[, c(depVariable)]
  #####################################################################
  variability <- list(
    "Value" = length(unique(columnOfInterest)) / max(length(columnOfInterest), 1),
    "Type" = "Total unique response divided by total number of response"
  )
  #####################################################################
  DSsize <- SummaryStats(
    x = x,
    formula = object$input$formula,
    # label = 'Summary statistics',
    lower = TRUE,
    drop = TRUE,
    sep = "_"
  )
  MultiBatch <- ifelse(multiBatch(x),
    "Dataset contains multiple batches",
    "Dataset contains single batch"
  )
  addInfo <- list(
    Data = list(
      "Data signature" = dataSignature(
        formula = frm,
        data = x
      ),
      "Variability" = variability,
      "Summary statistics" = DSsize
    ),
    Analysis = list(
      "Model setting" = extractFERRTerms(object),
      "Implementation specification" = RRextraDetailsExtractor(object),
      "Is model optimised" = NULL,
      "Multibatch in analysis" = MultiBatch,
      "Gender included in analysis" = GenderIncludedInAnalysis(x),
      "Further models" = ReFurtherModels(Vsplit),
      "Effect sizes" = "Look at the individual models",
      "Other residual normality tests" = NULL
    )
  )
  #####################################################################
  percentageChanges <- NULL
  #####################################################################
  OpenStatsReportRR0 <- list(
    "Applied method" = framework,
    "Dependent variable" = depVariable,
    "Batch included" = NULL,
    "Batch p-value" = NULL,
    "Residual variances homogeneity" = NULL,
    "Residual variances homogeneity p-value" = NULL,
    #####################################################################
    "Genotype contribution" = list(
      "Overall" = lowHighList(
        extractFisherSubTableResults1(VsplitLow$Genotype$result),
        extractFisherSubTableResults1(VsplitHig$Genotype$result),
        "Details"  = GenotypeDiscLabel
      ),
      "Sex FvKO p-value" = lowHighList(
        extractFisherSubTableResults1(VsplitLowFemale$Genotype$result),
        extractFisherSubTableResults1(VsplitHigFemale$Genotype$result),
        "Details" = SexDiscLabel
      ),
      "Sex MvKO p-value" = lowHighList(
        extractFisherSubTableResults1(VsplitLowMale$Genotype$result),
        extractFisherSubTableResults1(VsplitHigMale$Genotype$result),
        "Details" = SexDiscLabel
      ),
      "Sexual dimorphism detected" = list(
        "Criteria" = TRUE,
        "Note" = "Sex specific results for Low/High tables are always reported."
      )
    ),
    "Genotype estimate" = lowHighList(
      lapply1(VsplitLow$Genotype$result, CatEstimateAndCI),
      lapply1(VsplitHig$Genotype$result, CatEstimateAndCI),
      "Details" = GenotypeDiscLabel
    ),
    "Genotype standard error" = NULL,
    "Genotype p-value" = lowHighList(
      extractFisherSubTableResults1(VsplitLow$Genotype$result),
      extractFisherSubTableResults1(VsplitHig$Genotype$result),
      "Details" = GenotypeDiscLabel
    ),
    "Genotype percentage change" = percentageChanges,
    "Genotype effect size" = lowHighList(
      extractFisherSubTableResults1(VsplitLow$Genotype$result, "effect"),
      extractFisherSubTableResults1(VsplitHig$Genotype$result, "effect"),
      "Details" = GenotypeDiscLabel
    ),
    #####################################################################
    "Sex estimate" = lowHighList(
      lapply1(VsplitLow$Sex$result, CatEstimateAndCI),
      lapply1(VsplitHig$Sex$result, CatEstimateAndCI),
      "Details" = GenotypeDiscLabel
    ),
    "Sex standard error" = NULL,
    "Sex p-value" = lowHighList(
      extractFisherSubTableResults1(VsplitLow$Sex$result),
      extractFisherSubTableResults1(VsplitHig$Sex$result),
      "Details" = GenotypeDiscLabel
    ),
    "Sex effect size" = lowHighList(
      extractFisherSubTableResults1(VsplitLow$Sex$result, "effect"),
      extractFisherSubTableResults1(VsplitHig$Sex$result, "effect"),
      "Details" = GenotypeDiscLabel
    ),
    #####################################################################
    "LifeStage estimate" = lowHighList(
      lapply1(VsplitLow$LifeStage$result, CatEstimateAndCI),
      lapply1(VsplitHig$LifeStage$result, CatEstimateAndCI),
      "Details" = GenotypeDiscLabel
    ),
    "LifeStage standard error" = NULL,
    "LifeStage p-value" = lowHighList(
      extractFisherSubTableResults1(VsplitLow$LifeStage$result),
      extractFisherSubTableResults1(VsplitHig$LifeStage$result),
      "Details" = GenotypeDiscLabel
    ),
    "LifeStage effect size" = lowHighList(
      extractFisherSubTableResults1(VsplitLow$LifeStage$result, "effect"),
      extractFisherSubTableResults1(VsplitHig$LifeStage$result, "effect"),
      "Details" = GenotypeDiscLabel
    ),
    #####################################################################
    "Weight estimate" = NULL,
    "Weight standard error" = NULL,
    "Weight p-value" = NULL,
    "Weight effect size" = NULL,
    #####################################################################
    "Gp1 genotype" = Labels$Genotype$Control,
    "Gp1 Residuals normality test" = NULL,
    "Gp2 genotype" = Labels$Genotype$Mutant,
    "Gp2 Residuals normality test" = NULL,
    #####################################################################
    "Blups test" = NULL,
    "Rotated residuals normality test" = NULL,
    #####################################################################
    "Intercept estimate" = NULL,
    "Intercept standard error" = NULL,
    "Intercept p-value" = NULL,
    #####################################################################
    "Interactions included" = list(
      "Genotype Sex"                   =  NULL,
      "Genotype LifeStage"             =  NULL,
      "Sex LifeStage"                  =  NULL,
      "Genotype Sex LifeStage"         =  NULL
    ),
    #####################################################################
    ################ interaction
    "Interactions p-value" = list(
      "Genotype Sex"                  = NULL,
      "Genotype LifeStage"            = NULL,
      "Sex LifeStage"                 = NULL,
      "Genotype Sex LifeStage"        = NULL
    ),
    #####################################################################
    ################ Sex interactions
    "Sex FvKO estimate" = lowHighList(
      lapply1(VsplitLowFemale$Genotype$result, CatEstimateAndCI),
      lapply1(VsplitHigFemale$Genotype$result, CatEstimateAndCI),
      "Details" = SexDiscLabel
    ),
    "Sex FvKO standard error" = NULL,
    "Sex FvKO p-value" = lowHighList(
      extractFisherSubTableResults1(VsplitLowFemale$Genotype$result),
      extractFisherSubTableResults1(VsplitHigFemale$Genotype$result),
      "Details" = SexDiscLabel
    ),
    "Sex FvKO effect size" = lowHighList(
      extractFisherSubTableResults1(VsplitLowFemale$Genotype$result, "effect"),
      extractFisherSubTableResults1(VsplitHigFemale$Genotype$result, "effect"),
      "Details" = SexDiscLabel
    ),
    #####################################################################
    "Sex MvKO estimate" = lowHighList(
      lapply1(VsplitLowMale$Genotype$result, CatEstimateAndCI),
      lapply1(VsplitHigMale$Genotype$result, CatEstimateAndCI),
      "Details" = SexDiscLabel
    ),
    "Sex MvKO standard error" = NULL,
    "Sex MvKO p-value" = lowHighList(
      extractFisherSubTableResults1(VsplitLowMale$Genotype$result),
      extractFisherSubTableResults1(VsplitHigMale$Genotype$result),
      "Details" = SexDiscLabel
    ),
    "Sex MvKO effect size" = lowHighList(
      extractFisherSubTableResults1(VsplitLowMale$Genotype$result, "effect"),
      extractFisherSubTableResults1(VsplitHigMale$Genotype$result, "effect"),
      "Details" = SexDiscLabel
    ),
    #####################################################################
    ################ LifeStage interaction
    "LifeStage EvKO estimate" = lowHighList(
      lapply1(VsplitLowEarly$Genotype$result, CatEstimateAndCI),
      lapply1(VsplitHigEarly$Genotype$result, CatEstimateAndCI),
      "Details" = LifeStageDiscLabel
    ),
    "LifeStage EvKO standard error" = NULL,
    "LifeStage EvKO p-value" = lowHighList(
      extractFisherSubTableResults1(VsplitLowEarly$Genotype$result),
      extractFisherSubTableResults1(VsplitHigEarly$Genotype$result),
      "Details" = LifeStageDiscLabel
    ),
    "LifeStage EvKO effect size" = lowHighList(
      extractFisherSubTableResults1(VsplitLowEarly$Genotype$result, "effect"),
      extractFisherSubTableResults1(VsplitHigEarly$Genotype$result, "effect"),
      "Details" = LifeStageDiscLabel
    ),
    #####################################################################
    "LifeStage LvKO estimate" = lowHighList(
      lapply1(VsplitLowLate$Genotype$result, CatEstimateAndCI),
      lapply1(VsplitHigLate$Genotype$result, CatEstimateAndCI),
      "Details" = LifeStageDiscLabel
    ),
    "LifeStage LvKO standard error" = NULL,
    "LifeStage LvKO p-value" = lowHighList(
      extractFisherSubTableResults1(VsplitLowLate$Genotype$result),
      extractFisherSubTableResults1(VsplitHigLate$Genotype$result),
      "Details" = LifeStageDiscLabel
    ),
    "LifeStage LvKO effect size" = lowHighList(
      extractFisherSubTableResults1(VsplitLowLate$Genotype$result, "effect"),
      extractFisherSubTableResults1(VsplitHigLate$Genotype$result, "effect"),
      "Details" = LifeStageDiscLabel
    ),
    #####################################################################
    ################ Sex LifeStage Genotype interactions
    #####################################################################
    # 1.
    "LifeStageSexGenotype FvEvKO estimate" = lowHighList(
      lapply1(VsplitLowEarlyFemale$Genotype$result, CatEstimateAndCI),
      lapply1(VsplitHigEarlyFemale$Genotype$result, CatEstimateAndCI),
      "Details" = LifeStageSexDiscLabel
    ),
    "LifeStageSexGenotype FvEvKO standard error" = NULL,
    "LifeStageSexGenotype FvEvKO p-value" = lowHighList(
      extractFisherSubTableResults1(VsplitLowEarlyFemale$Genotype$result),
      extractFisherSubTableResults1(VsplitHigEarlyFemale$Genotype$result),
      "Details" = LifeStageSexDiscLabel
    ),
    "LifeStageSexGenotype FvEvKO effect size" = lowHighList(
      extractFisherSubTableResults1(VsplitLowEarlyFemale$Genotype$result, "effect"),
      extractFisherSubTableResults1(VsplitHigEarlyFemale$Genotype$result, "effect"),
      "Details" = LifeStageSexDiscLabel
    ),
    # 2.
    "LifeStageSexGenotype MvEvKO estimate" = lowHighList(
      lapply1(VsplitLowEarlyMale$Genotype$result, CatEstimateAndCI),
      lapply1(VsplitHigEarlyMale$Genotype$result, CatEstimateAndCI),
      "Details" = LifeStageSexDiscLabel
    ),
    "LifeStageSexGenotype MvEvKO standard error" = NULL,
    "LifeStageSexGenotype MvEvKO p-value" = lowHighList(
      extractFisherSubTableResults1(VsplitLowEarlyMale$Genotype$result),
      extractFisherSubTableResults1(VsplitHigEarlyMale$Genotype$result),
      "Details" = LifeStageSexDiscLabel
    ),
    "LifeStageSexGenotype MvEvKO effect size" = lowHighList(
      extractFisherSubTableResults1(VsplitLowEarlyMale$Genotype$result, "effect"),
      extractFisherSubTableResults1(VsplitHigEarlyMale$Genotype$result, "effect"),
      "Details" = LifeStageSexDiscLabel
    ),
    # 3.
    "LifeStageSexGenotype FvLvKO estimate" = lowHighList(
      lapply1(VsplitLowLateFemale$Genotype$result, CatEstimateAndCI),
      lapply1(VsplitHigLateFemale$Genotype$result, CatEstimateAndCI),
      "Details" = LifeStageSexDiscLabel
    ),
    "LifeStageSexGenotype FvLvKO standard error" = NULL,
    "LifeStageSexGenotype FvLvKO p-value" = lowHighList(
      extractFisherSubTableResults1(VsplitLowLateFemale$Genotype$result),
      extractFisherSubTableResults1(VsplitHigLateFemale$Genotype$result),
      "Details" = LifeStageSexDiscLabel
    ),
    "LifeStageSexGenotype FvLvKO effect size" = lowHighList(
      extractFisherSubTableResults1(VsplitLowLateFemale$Genotype$result, "effect"),
      extractFisherSubTableResults1(VsplitHigLateFemale$Genotype$result, "effect"),
      "Details" = LifeStageSexDiscLabel
    ),
    "LifeStageSexGenotype MvLvKO estimate" = lowHighList(
      lapply1(VsplitLowLateMale$Genotype$result, CatEstimateAndCI),
      lapply1(VsplitHigLateMale$Genotype$result, CatEstimateAndCI),
      "Details" = LifeStageSexDiscLabel
    ),
    "LifeStageSexGenotype MvLvKO standard error" = NULL,
    "LifeStageSexGenotype MvLvKO p-value" = lowHighList(
      extractFisherSubTableResults1(VsplitLowLateMale$Genotype$result),
      extractFisherSubTableResults1(VsplitHigLateMale$Genotype$result),
      "Details" = LifeStageSexDiscLabel
    ),
    "LifeStageSexGenotype MvLvKO effect size" = lowHighList(
      extractFisherSubTableResults1(VsplitLowLateMale$Genotype$result, "effect"),
      extractFisherSubTableResults1(VsplitHigLateMale$Genotype$result, "effect"),
      "Details" = LifeStageSexDiscLabel
    ),
    ################
    "Classification tag" = NULL,
    "Transformation" = NULL,
    "Additional information" = addInfo
  )

  return(OpenStatsReportRR0)
}
