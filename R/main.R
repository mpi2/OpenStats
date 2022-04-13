OpenStatsAnalysis <- function(OpenStatsListObject = NULL,
                              method = NULL,
                              MM_fixed = TypicalModel(
                                depVariable = "data_point",
                                withWeight = MM_BodyWeightIncluded,
                                Sex = TRUE,
                                LifeStage = TRUE,
                                data = OpenStatsListObject@datasetPL,
                                others = NULL,
                                debug = debug
                              ),
                              MM_random = rndProce("TYPICAL"),
                              MM_BodyWeightIncluded = TRUE,
                              MM_lower = ~ Genotype + 1,
                              MM_weight = if (TermInModelAndnLevels(
                                model = MM_fixed,
                                data = OpenStatsListObject@datasetPL
                              )) {
                                varIdent(form = ~ 1 |
                                  LifeStage)
                              } else {
                                varIdent(form = ~ 1 | Genotype)
                              },
                              MM_direction = "both",
                              MM_checks = c(TRUE, TRUE, TRUE, TRUE),
                              MM_optimise = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
                              ### FE or RR
                              FE_formula = category ~ Genotype + Sex + LifeStage,
                              RR_formula = data_point ~ Genotype + Sex + LifeStage,
                              RRrefLevel = "control",
                              RR_prop = 0.95,
                              FERR_rep = 1500,
                              FERR_FullComparisions = c(TRUE, FALSE),
                              ##### Others
                              MMFERR_conf.level = 0.95,
                              debug = TRUE,
                              ...) {
  r <- NULL
  s <- tryCatch(
    expr = {
      suppressMessagesANDWarnings(
        OpenStatsAnalysis0(
          OpenStatsListObject = OpenStatsListObject,
          method = method,
          MM_fixed = MM_fixed,
          MM_random = MM_random,
          MM_lower = MM_lower,
          MM_weight = MM_weight,
          MM_direction = MM_direction,
          MM_checks   = MM_checks,
          MM_optimise = MM_optimise,
          FE_formula = FE_formula,
          RR_formula = RR_formula,
          RR_prop = RR_prop,
          RRrefLevel = RRrefLevel,
          FERR_rep = FERR_rep,
          FERR_FullComparisions = FERR_FullComparisions,
          debug = debug,
          MMFERR_conf.level = MMFERR_conf.level,
          ...
        ),
        sup.messages = !debug,
        sup.warnings = TRUE
      )
    },
    warning = function(war) {
      message0("This operation failed with a warning (see below): ")
      r$messages$warning <- war
      warning(war)
      return(NULL)
    },
    error = function(err) {
      message0("This operation failed with an error (see below): ")
      r$messages$error <- err
      message(err)
      return(NULL)
    }
  )
  if (is.null(s)) {
    return(r)
  } else {
    return(s)
  }
}

OpenStatsAnalysis0 <- function(OpenStatsListObject = NULL,
                               method,
                               MM_fixed,
                               MM_random,
                               MM_lower,
                               MM_weight,
                               MM_direction = "both",
                               MM_checks,
                               MM_optimise,
                               FE_formula,
                               RR_formula,
                               RR_prop,
                               RRrefLevel,
                               FERR_rep,
                               FERR_FullComparisions,
                               MMFERR_conf.level = 0.95,
                               debug = TRUE,
                               ...) {
  message0("OpenStats loaded.")
  if (is.null(OpenStatsListObject)) {
    stop("\n ~> The input dataset cannot be NULL")
  }
  if (!is0(OpenStatsListObject, "PhenList") &&
    !is0(OpenStatsListObject, "OpenStatsList")) {
    stop('\n ~> function expects "PhenList" or "OpenStatsList" object \n')
  }
  if (noVariation(data = OpenStatsListObject@datasetPL)) {
    stop("\n ~> There is no variation in Genotype.\n")
  }
  if (MMFERR_conf.level >= 1 ||
    MMFERR_conf.level <= 0 ||
    length(MMFERR_conf.level) > 1) {
    stop("\n ~> Confidence level must be a single value in (0,1) interval")
  }
  if (length(method) > 1 || !method %in% c("MM", "FE", "RR")) {
    stop("\n ~> `method` must be one of `MM`, `FE` or `RR`")
  }
  #####
  if (toupper(method) %in% "MM") {
    if (length(MM_optimise) != 6) {
      stop("\n ~> `MM_optimise` must be a vector of 6 TRUE/FALSE elements. Example:\n\t c(TRUE,TRUE,TRUE,TRUE,TRUE,TRUE)")
    }
    if (FormulaContainsFunction(MM_fixed)) {
      stop(
        "\n ~> Function detected in `MM_fixed`.\n\t	Please apply the function to the dataset prior to the analysis."
      )
    }
    message0("Linear Mixed Model (MM framework) in progress ...")
    output <- M.opt(
      fixed     = MM_fixed,
      random    = MM_random,
      object    = OpenStatsListObject,
      lower     = MM_lower,
      direction = MM_direction,
      weight    = MM_weight,
      checks    = MM_checks,
      optimise  = MM_optimise,
      trace     = FALSE,
      method    = "MM",
      ci_levels = MMFERR_conf.level,
      ...
    )
    # Important!
    if (!is.null(output$input)) {
      output$input$fixed <- MM_fixed
    }
  } else if (toupper(method) %in% "FE") {
    if (FormulaContainsFunction(FE_formula)) {
      stop(
        "\n ~> Function detected in `FE_formula`.\n\t	Please apply the function to the dataset prior to the analysis."
      )
    }
    if (length(FERR_FullComparisions) != 2) {
      stop(
        "\n ~> `FERR_FullComparisions` must be a vector of 2 TRUE/FALSE elements. Example:\n\t c(TRUE,TRUE)"
      )
    }
    message0("Fisher Exact Test (FE framework) in progress ...")
    FERR_FullComparisionsMessage(FERR_FullComparisions)
    output <- crunner(
      object = OpenStatsListObject,
      formula = MoveResponseToRightOfTheFormula(FE_formula),
      # expandDottedFormula(formula = FE_formula, data = OpenStatsListObject@datasetPL)  ,
      rep              = FERR_rep,
      method           = "FE",
      fullComparisions = FERR_FullComparisions[1],
      ci_levels        = MMFERR_conf.level,
      InterLevelComparisions = FERR_FullComparisions[2],

      ...
    )
    # Important!
    if (!is.null(output$input)) {
      output$input$formula <- FE_formula
    }
  } else if (toupper(method) %in% "RR") {
    if (FormulaContainsFunction(RR_formula)) {
      stop(
        "\n ~> Function detected in `RR_formula`.\n\t	Please apply the function to the dataset prior to the analysis."
      )
    }
    if (length(FERR_FullComparisions) != 2) {
      stop(
        "\n ~> `FERR_FullComparisions` must be a vector of 2 TRUE/FALSE elements. Example:\n\t c(TRUE,TRUE)"
      )
    }
    message0("Reference Range Plus (RR framework) in progress ...")
    FERR_FullComparisionsMessage(FERR_FullComparisions)
    output <- RRrunner(
      object  = OpenStatsListObject,
      formula = MoveResponseToRightOfTheFormula(RR_formula),
      rep = FERR_rep,
      method = "RR",
      RRprop = RR_prop,
      ci_levels = MMFERR_conf.level,
      RRrefLevel = RRrefLevel,
      fullComparisions = FERR_FullComparisions[1],
      InterLevelComparisions = FERR_FullComparisions[2],
      ...
    )
    # Important!
    if (!is.null(output$input)) {
      output$input$formula <- RR_formula
    }
  } else {
    message0('No "method" is specified. ')
    output <- NULL
  }
  return(output)
}
