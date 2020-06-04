OpenStatsReportCat <- function(object) {
  if (!is.null(object$messages)) {
    return(NULL)
  }
  #####################################################################
  Labels <- OpenStatsListLevels(object = object)
  Fmodel <- object$extra$Cleanedformula
  frm <- formula(Fmodel)
  depVariable <- all_vars0(frm)[1]
  # 	equation       = NULL
  formula <- printformula(frm)
  framework <- "Fisher Exact Test framework"
  # 	fittingMethod  = NULL
  #####################################################################
  x <- object$input$data
  columnOfInterest <- x[, c(depVariable)]
  objRes <- object$output$SplitModels[[depVariable]]
  #####################################################################
  variability <- list(
    "Value" = length(unique(columnOfInterest)) / max(length(columnOfInterest), 1),
    "Type" = "Total unique response divided by total number of response"
  )
  #####################################################################
  DSsize <- SummaryStats(
    x = x,
    formula = object$input$formula,
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
      # 'Formula'                = list(
      # 	input   = printformula(object$input$formula),
      # 	final   = printformula(formula)
      # ),
      "Model setting" = extractFERRTerms(object),
      "Is model optimised" = NULL,
      "Multibatch in analysis" = MultiBatch,
      "Gender included in analysis" = GenderIncludedInAnalysis(x),
      "Further models" = FeFurtherModels(object), 
      "Effect sizes" = "Look at the individual models",
      "Other residual normality tests" = NULL
    )
  )
  #####################################################################
  percentageChanges <- NULL
  #####################################################################
  OpenStatsReportFE0 <- list(
    "Applied method" = framework,
    "Dependent variable" = depVariable,
    "Batch included" = NULL,
    "Batch p-value" = NULL,
    "Residual variances homogeneity" = NULL,
    "Residual variances homogeneity p-value" = NULL,
    #####################################################################
    "Genotype contribution" = list(
      "Overall" = extractFisherSubTableResults(objRes$Genotype$result),
      "Sex FvKO p-value" = extractFisherSubTableResults(objRes$Genotype_Female$result),
      "Sex MvKO p-value" = extractFisherSubTableResults(objRes$Genotype_Male$result),
      "Sexual dimorphism detected" = list(
        "Criteria" = TRUE,
        "Note" = "Sex specific results are always reported."
      )
    ),
    "Genotype estimate" = lapply0(objRes$Genotype$result, CatEstimateAndCI),
    "Genotype standard error" = NULL,
    "Genotype p-value" = extractFisherSubTableResults(objRes$Genotype$result),
    "Genotype percentage change" = percentageChanges,
    "Genotype effect size" = extractFisherSubTableResults(objRes$Genotype$result, "effect"),
    #####################################################################
    "Sex estimate" = lapply0(objRes$Sex$result, CatEstimateAndCI),
    "Sex standard error" = NULL,
    "Sex p-value" = extractFisherSubTableResults(objRes$Sex$result),
    "Sex effect size" = extractFisherSubTableResults(objRes$Sex$result, "effect"),
    #####################################################################
    "LifeStage estimate" = lapply0(objRes$LifeStage$result, CatEstimateAndCI),
    "LifeStage standard error" = NULL,
    "LifeStage p-value" = extractFisherSubTableResults(objRes$LifeStage$result),
    "LifeStage effect size" = extractFisherSubTableResults(objRes$LifeStage$result, "effect"),
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
      "Genotype Sex"            = NULL,
      "Genotype LifeStage"      = NULL,
      "Sex LifeStage"           = NULL,
      "Genotype Sex LifeStage"  = NULL
    ),
    #####################################################################
    ################ Sex interactions
    "Sex FvKO estimate" = lapply0(objRes$Genotype_Female$result, CatEstimateAndCI),
    "Sex FvKO standard error" = NULL,
    "Sex FvKO p-value" = extractFisherSubTableResults(objRes$Genotype_Female$result),
    "Sex FvKO effect size" = extractFisherSubTableResults(objRes$Genotype_Female$result, "effect"),
    #####################################################################
    "Sex MvKO estimate" = lapply0(objRes$Genotype_Male$result, CatEstimateAndCI),
    "Sex MvKO standard error" = NULL,
    "Sex MvKO p-value" = extractFisherSubTableResults(objRes$Genotype_Male$result),
    "Sex MvKO effect size" = extractFisherSubTableResults(objRes$Genotype_Male$result, "effect"),
    #####################################################################
    ################ LifeStage interaction
    "LifeStage EvKO estimate" = lapply0(objRes$Genotype_Early$result, CatEstimateAndCI),
    "LifeStage EvKO standard error" = NULL,
    "LifeStage EvKO p-value" = extractFisherSubTableResults(objRes$Genotype_Early$result),
    "LifeStage EvKO effect size" = extractFisherSubTableResults(objRes$Genotype_Early$result, "effect"),
    #####################################################################
    "LifeStage LvKO estimate" = lapply0(objRes$Genotype_Late$result, CatEstimateAndCI),
    "LifeStage LvKO standard error" = NULL,
    "LifeStage LvKO p-value" = extractFisherSubTableResults(objRes$Genotype_Late$result),
    "LifeStage LvKO effect size" = extractFisherSubTableResults(objRes$Genotype_Late$result, "effect"),
    #####################################################################
    ################ Sex LifeStage Genotype interactions
    # 1.
    "LifeStageSexGenotype FvEvKO estimate" = lapply0(objRes$Genotype_Female.Early$result, CatEstimateAndCI),
    "LifeStageSexGenotype FvEvKO standard error" = NULL,
    "LifeStageSexGenotype FvEvKO p-value" = extractFisherSubTableResults(objRes$Genotype_Female.Early$result),
    "LifeStageSexGenotype FvEvKO effect size" = extractFisherSubTableResults(objRes$Genotype_Female.Early$result, "effect"),
    # 2.
    "LifeStageSexGenotype MvEvKO estimate" = lapply0(objRes$Genotype_Male.Early$result, CatEstimateAndCI),
    "LifeStageSexGenotype MvEvKO standard error" = NULL,
    "LifeStageSexGenotype MvEvKO p-value" = extractFisherSubTableResults(objRes$Genotype_Male.Early$result),
    "LifeStageSexGenotype MvEvKO effect size" = extractFisherSubTableResults(objRes$Genotype_Male.Early$result, "effect"),
    # 3.
    "LifeStageSexGenotype FvLvKO estimate" = lapply0(objRes$Genotype_Female.Late$result, CatEstimateAndCI),
    "LifeStageSexGenotype FvLvKO standard error" = NULL,
    "LifeStageSexGenotype FvLvKO p-value" = extractFisherSubTableResults(objRes$Genotype_Female.Late$result),
    "LifeStageSexGenotype FvLvKO effect size" = extractFisherSubTableResults(objRes$Genotype_Female.Late$result, "effect"),

    "LifeStageSexGenotype MvLvKO estimate" = lapply0(objRes$Genotype_Male.Late$result, CatEstimateAndCI),
    "LifeStageSexGenotype MvLvKO standard error" = NULL,
    "LifeStageSexGenotype MvLvKO p-value" = extractFisherSubTableResults(objRes$Genotype_Male.Late$result),
    "LifeStageSexGenotype MvLvKO effect size" = extractFisherSubTableResults(objRes$Genotype_Male.Late$result, "effect"),
    ################
    "Classification tag" = NULL,
    "Transformation" = NULL,
    "Additional information" = addInfo
  )
  return(OpenStatsReportFE0)
}
