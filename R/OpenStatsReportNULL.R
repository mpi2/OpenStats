OpenStatsReportNULL <- function(object) {
  VO <- list(
    "Applied method" = NULL,
    "Dependent variable" = NULL,
    "Batch included" = NULL,
    "Batch p-value" = NULL,
    "Residual variances homogeneity" = NULL,
    "Residual variances homogeneity p-value" = NULL,
    #####################################################################
    "Genotype contribution" = list(
      "Overall"            = NULL,
      "Sex FvKO p-value"   = NULL,
      "Sex MvKO p-value"   = NULL,
      "Sexual dimorphism detected" = NULL
    ),
    "Genotype estimate" = NULL,
    "Genotype standard error" = NULL,
    "Genotype p-value" = NULL,
    "Genotype percentage change" = NULL,
    "Genotype effect size" = NULL,
    #####################################################################
    "Sex estimate" = NULL,
    "Sex standard error" = NULL,
    "Sex p-value" = NULL,
    "Sex effect size" = NULL,
    #####################################################################
    "LifeStage estimate" = NULL,
    "LifeStage standard error" = NULL,
    "LifeStage p-value" = NULL,
    "LifeStage effect size" = NULL,
    #####################################################################
    "Weight estimate" = NULL,
    "Weight standard error" = NULL,
    "Weight p-value" = NULL,
    "Weight effect size" = NULL,
    #####################################################################
    "Gp1 genotype" = NULL,
    "Gp1 Residuals normality test" = NULL,
    "Gp2 genotype" = NULL,
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
    "Sex FvKO estimate" = NULL,
    "Sex FvKO standard error" = NULL,
    "Sex FvKO p-value" = NULL,
    "Sex FvKO effect size" = NULL,
    #####################################################################
    "Sex MvKO estimate" = NULL,
    "Sex MvKO standard error" = NULL,
    "Sex MvKO p-value" = NULL,
    "Sex MvKO effect size" = NULL,
    #####################################################################
    ################ LifeStage interaction
    "LifeStage EvKO estimate" = NULL,
    "LifeStage EvKO standard error" = NULL,
    "LifeStage EvKO p-value" = NULL,
    "LifeStage EvKO effect size" = NULL,
    #####################################################################
    "LifeStage LvKO estimate" = NULL,
    "LifeStage LvKO standard error" = NULL,
    "LifeStage LvKO p-value" = NULL,
    "LifeStage LvKO effect size" = NULL,
    #####################################################################
    ################ Sex LifeStage Genotype interactions
    # 1.
    "LifeStageSexGenotype FvEvKO estimate" = NULL,
    "LifeStageSexGenotype FvEvKO standard error" = NULL,
    "LifeStageSexGenotype FvEvKO p-value" = NULL,
    "LifeStageSexGenotype FvEvKO effect size" = NULL,
    # 2.
    "LifeStageSexGenotype MvEvKO estimate" = NULL,
    "LifeStageSexGenotype MvEvKO standard error" = NULL,
    "LifeStageSexGenotype MvEvKO p-value" = NULL,
    "LifeStageSexGenotype MvEvKO effect size" = NULL,
    # 3.
    "LifeStageSexGenotype FvLvKO estimate" = NULL,
    "LifeStageSexGenotype FvLvKO standard error" = NULL,
    "LifeStageSexGenotype FvLvKO p-value" = NULL,
    "LifeStageSexGenotype FvLvKO effect size" = NULL,

    "LifeStageSexGenotype MvLvKO estimate" = NULL,
    "LifeStageSexGenotype MvLvKO standard error" = NULL,
    "LifeStageSexGenotype MvLvKO p-value" = NULL,
    "LifeStageSexGenotype MvLvKO effect size" = NULL,
    ################
    "Classification tag" = NULL,
    "Transformation" = NULL,
    "Additional information" = object$messages
  )
  return(VO)
}
