classificationTag <- function(object = NULL,
                              phenotypeThreshold = 0.0001,
                              SexSpecificThreshold = .05,
                              userMode = "summaryOutput") {
  decimals <- ifelse(phenotypeThreshold < 1, 1, 0) * abs(log(phenotypeThreshold, 10))
  if (is.null(object) ||
    !is.null(object$messages)) {
    message0("~> The input object is missing or a failure object.")
    return(NULL)
  }

  if (!(userMode %in% c("summaryOutput", "vectorOutput"))) {
    stop(
      "~> Please define `userMode` you would like to use from `summaryOutout` or `vectorOutput`."
    )
  }
  if (length(phenotypeThreshold) != 1 ||
    phenotypeThreshold > 1 ||
    phenotypeThreshold < 0) {
    stop("~> `phenotypeThreshold` must be a single value in [0,1] interval.")
  }
  ######################
  debug <- FALSE
  ChangeClassification <- NA
  OpenStatListObj <- object$input$OpenStatsList
  tbl <- apply(
    table(
      OpenStatListObj@datasetUNF[, OpenStatListObj@dataset.colname.sex],
      OpenStatListObj@datasetUNF[, OpenStatListObj@dataset.colname.genotype]
    ),
    1,
    prod
  )
  csex <- asFactorAndSelectVariable(
    OpenStatListObj@datasetUNF,
    OpenStatListObj@dataset.colname.sex
  )
  SexInTheCleanedInputModel <- "Sex" %in% all_vars0(object$extra$Cleanedformula)
  if (!SexInTheCleanedInputModel) {
    csex <- NULL
  }
  nsex <- sum(tbl > 0) * !is.null(csex) # nlevels(csex)
  lsex <- names(tbl)[tbl > 0] # pasteComma(levels(csex))
  v <- OpenStatsReport(object)
  Labels <- OpenStatsListLevels(object = object)
  SexInteraction <- NULL
  ######################
  if (object$input$method %in% c("MM")) {
    SexInteraction <- TermInFormulaReturn(
      active = TRUE,
      formula = formula(object$output$Final.Model),
      term = CombineLevels(Labels$Sex$Sex,
        Labels$Genotype$Genotype,
        debug = debug
      ),
      not = NULL,
      return = modelSummaryPvalueExtract(
        x = object$output$SplitModels$Genotype_Sex,
        variable = CombineLevels(Labels$Sex$Sex,
          Labels$Genotype$Genotype,
          debug = debug
        ),
        anova = TRUE,
        debug = debug
      ),
      debug = debug
    )
    ##########################
    if (userMode == "summaryOutput") {
      if (length(v$`Genotype p-value`) < 1) {
        ChangeClassification <- NA
      } else if (v$`Genotype p-value` > phenotypeThreshold &&
        ifelse(length(SexInteraction) > 0,
          SexInteraction > phenotypeThreshold,
          TRUE
        )) {
        if (nsex == 1) {
          ChangeClassification <- paste0(
            "With phenotype threshold value ",
            phenotypeThreshold,
            " - no significant change for the one sex (",
            lsex,
            ") tested"
          )
        } else {
          ChangeClassification <- paste(
            "With phenotype threshold value",
            phenotypeThreshold,
            "- no significant change"
          )
        }
      } else {
        if (length(SexInteraction) < 1) {
          if (nsex == 1) {
            ChangeClassification <-
              paste0(
                "With phenotype threshold value ",
                phenotypeThreshold,
                " - a significant change for the one sex (",
                lsex,
                ") tested"
              )
          } else if (nsex == 2) {
            ChangeClassification <- paste(
              "With phenotype threshold value",
              phenotypeThreshold,
              "- both sexes equally"
            )
          } else {
            ChangeClassification <- paste(
              "With phenotype threshold value",
              phenotypeThreshold,
              "- regardless of gender"
            )
          }
        } else if (v$`Sex FvKO p-value` >= SexSpecificThreshold &&
          v$`Sex MvKO p-value` >= SexSpecificThreshold) {
          ChangeClassification <- paste0(
            "With phenotype threshold value ",
            phenotypeThreshold,
            " - cannot classify effect [",
            "Interaction pvalue = ",
            FormatPvalueClassificationTag(SexInteraction, decimals = decimals),
            ", Genotype Female pvalue = ",
            FormatPvalueClassificationTag(v$`Sex FvKO p-value`, decimals = decimals),
            ", Genotype Male pvalue = ",
            FormatPvalueClassificationTag(v$`Sex MvKO p-value`, decimals = decimals),
            "]"
          )
        } else if (v$`Sex FvKO p-value` < SexSpecificThreshold &&
          v$`Sex MvKO p-value` >= SexSpecificThreshold) {
          ChangeClassification <- paste(
            "With phenotype threshold value",
            phenotypeThreshold,
            "- females only"
          )
        } else if (v$`Sex FvKO p-value` >= SexSpecificThreshold &&
          v$`Sex MvKO p-value` < SexSpecificThreshold) {
          ChangeClassification <- paste(
            "With phenotype threshold value",
            phenotypeThreshold,
            "- males only"
          )
        } else if (as.list(v$`Sex FvKO estimate`)$Value > 0 &&
          as.list(v$`Sex MvKO estimate`)$Value > 0 ||
          as.list(v$`Sex FvKO estimate`)$Value < 0 &&
            as.list(v$`Sex MvKO estimate`)$Value < 0) {
          if (abs(as.list(v$`Sex FvKO estimate`)$Value) > abs(as.list(v$`Sex MvKO estimate`)$Value)) {
            ChangeClassification <-
              paste(
                "With phenotype threshold value",
                phenotypeThreshold,
                "- different size as females greater"
              )
          } else {
            ChangeClassification <- paste(
              "With phenotype threshold value",
              phenotypeThreshold,
              "- different size as males greater"
            )
          }
        } else {
          ChangeClassification <- paste(
            "With phenotype threshold value",
            phenotypeThreshold,
            "- different direction for the sexes"
          )
        }
      }
    } else {
      if (length(SexInteraction) < 1) {
        if (nsex == 1) {
          ChangeClassification <-
            paste0(
              "If phenotype is significant it is for the one sex (",
              lsex,
              ") tested"
            )
        } else if (nsex == 2) {
          ChangeClassification <-
            paste("If phenotype is significant - both sexes equally")
        } else {
          ChangeClassification <-
            paste("If phenotype is significant - regardless of gender")
        }
      } else
      if (v$`Sex FvKO p-value` >= SexSpecificThreshold &&
        v$`Sex MvKO p-value` >= SexSpecificThreshold) {
        ChangeClassification <-
          paste0(
            "If phenotype is significant ",
            "- cannot classify effect [",
            "Interaction pvalue = ",
            FormatPvalueClassificationTag(SexInteraction, decimals = decimals),
            ", Genotype Female pvalue = ",
            FormatPvalueClassificationTag(v$`Sex FvKO p-value`, decimals = decimals),
            ", Genotype Male pvalue = ",
            FormatPvalueClassificationTag(v$`Sex MvKO p-value`, decimals = decimals),
            "]"
          )
      } else if (v$`Sex FvKO p-value` < SexSpecificThreshold &&
        v$`Sex MvKO p-value` >= SexSpecificThreshold) {
        ChangeClassification <-
          paste("If phenotype is significant - females only")
      } else if (v$`Sex FvKO p-value` >= SexSpecificThreshold &&
        v$`Sex MvKO p-value` < SexSpecificThreshold) {
        ChangeClassification <-
          paste("If phenotype is significant - males only")
      } else
      if (as.list(v$`Sex FvKO estimate`)$Value > 0 &&
        as.list(v$`Sex MvKO estimate`)$Value > 0 ||
        as.list(v$`Sex FvKO estimate`)$Value < 0 &&
          as.list(v$`Sex MvKO estimate`)$Value < 0) {
        if (abs(as.list(v$`Sex FvKO estimate`)$Value) > abs(as.list(v$`Sex MvKO estimate`)$Value)) {
          ChangeClassification <-
            paste("If phenotype is significant - different size as females greater")
        } else {
          ChangeClassification <-
            paste("If phenotype is significant - different size as males greater")
        }
      } else {
        ChangeClassification <-
          paste("If phenotype is significant - different direction for the sexes")
      }
    }
  } else if (object$input$method == "FE") {
    if (userMode == "summaryOutput") {
      all_p.value <- NullOrValue(v$`Genotype contribution`$`Overall`$`Complete table`$p.value)
      female_p.value <- NullOrValue(v$`Genotype contribution`$`Sex FvKO p-value`$`Complete table`$p.value)
      male_p.value <- NullOrValue(v$`Genotype contribution`$`Sex MvKO p-value`$`Complete table`$p.value)

      if (nsex == 1) {
        ChangeClassification <-
          paste0("Not significant for the one sex (", lsex, ") tested")
      } else {
        ChangeClassification <- paste("Not significant")
      }

      # Tag
      # combined & males & females
      if (all_p.value < phenotypeThreshold &&
        male_p.value < phenotypeThreshold &&
        female_p.value < phenotypeThreshold) {
        ChangeClassification <-
          paste(
            "With phenotype threshold value",
            phenotypeThreshold,
            "- significant in males, females and in combined dataset"
          )
      }
      # combined & males & !females
      if (all_p.value < phenotypeThreshold &&
        male_p.value < phenotypeThreshold &&
        female_p.value >= phenotypeThreshold) {
        ChangeClassification <-
          paste(
            "With phenotype threshold value",
            phenotypeThreshold,
            "- significant in males and in combined dataset"
          )
      }
      # combined & !males & females
      if (all_p.value < phenotypeThreshold &&
        male_p.value >= phenotypeThreshold &&
        female_p.value < phenotypeThreshold) {
        ChangeClassification <-
          paste(
            "With phenotype threshold value",
            phenotypeThreshold,
            "- significant in females and in combined dataset"
          )
      }
      # combined & !males & !females
      if (all_p.value < phenotypeThreshold &&
        male_p.value >= phenotypeThreshold &&
        female_p.value >= phenotypeThreshold) {
        if (nsex == 2) {
          ChangeClassification <-
            paste(
              "With phenotype threshold value",
              phenotypeThreshold,
              "- significant in combined dataset only"
            )
        } else {
          ChangeClassification <-
            paste0(
              "With phenotype threshold value ",
              phenotypeThreshold,
              " - significant for the sex (",
              lsex,
              ") tested"
            )
        }
      }
      # !combined & males & females
      if (all_p.value >= phenotypeThreshold &&
        male_p.value < phenotypeThreshold &&
        female_p.value < phenotypeThreshold) {
        ChangeClassification <-
          paste(
            "With phenotype threshold value",
            phenotypeThreshold,
            "- significant in males and in females datasets"
          )
      }
      # !combined & males & !females
      if (all_p.value >= phenotypeThreshold &&
        male_p.value < phenotypeThreshold &&
        female_p.value >= phenotypeThreshold) {
        ChangeClassification <-
          paste(
            "With phenotype threshold value",
            phenotypeThreshold,
            "- significant in males dataset only"
          )
      }
      # !combined & !males & females
      if (all_p.value >= phenotypeThreshold &&
        male_p.value >= phenotypeThreshold &&
        female_p.value < phenotypeThreshold) {
        ChangeClassification <-
          paste(
            "With phenotype threshold value",
            phenotypeThreshold,
            "- significant in females dataset only"
          )
      }
    }
    else {
      ChangeClassification <- NA
    }
  } else if (object$input$method == "RR") {
    if (userMode == "summaryOutput") {
      direction_all <- NA
      all_p.value <- 10
      direction_females <- NA
      female_p.value <- 10
      direction_males <- NA
      male_p.value <- 10
      high_male_p.value <- NullOrValue(v$`Genotype contribution`$`Sex MvKO p-value`$High$p.value)
      high_female_p.value <- NullOrValue(v$`Genotype contribution`$`Sex FvKO p-value`$High$p.value)
      high_all_p.value <- NullOrValue(v$`Genotype contribution`$`Overall`$High$p.value)
      low_male_p.value <- NullOrValue(v$`Genotype contribution`$`Sex MvKO p-value`$Low$p.value)
      low_female_p.value <- NullOrValue(v$`Genotype contribution`$`Sex FvKO p-value`$Low$p.value)
      low_all_p.value <- NullOrValue(v$`Genotype contribution`$Overall$Low$p.value)


      # High classification p-val is less than threshold and low classification p-val is more that threshold
      if (high_all_p.value < phenotypeThreshold &&
        low_all_p.value >= phenotypeThreshold) {
        direction_all <- "High"
        all_p.value <- high_all_p.value
      } else if (high_all_p.value >= phenotypeThreshold &&
        low_all_p.value < phenotypeThreshold) {
        # Low classification p-val is less than threshold and high classification p-val is more that threshold
        direction_all <- "Low"
        all_p.value <- low_all_p.value
      }

      if (high_male_p.value < phenotypeThreshold &&
        low_male_p.value >= phenotypeThreshold) {
        direction_males <- "High"
        male_p.value <- high_male_p.value
      }
      else if (high_male_p.value >= phenotypeThreshold &&
        low_male_p.value < phenotypeThreshold) {
        direction_males <- "Low"
        male_p.value <- low_male_p.value
      }


      if (high_female_p.value < phenotypeThreshold &&
        low_female_p.value >= phenotypeThreshold) {
        direction_females <- "High"
        female_p.value <- high_female_p.value
      }
      else if (high_female_p.value >= phenotypeThreshold &&
        low_female_p.value < phenotypeThreshold) {
        direction_females <- "Low"
        female_p.value <- low_female_p.value
      }

      if (nsex == 1) {
        ChangeClassification <-
          paste0("Not significant for the one sex (", lsex, ") tested")
      }
      else {
        ChangeClassification <- paste("Not significant")
      }

      if (nsex == 2) {
        # Tag
        # combined & males & females
        if (all_p.value < phenotypeThreshold &&
          male_p.value < phenotypeThreshold
        && female_p.value < phenotypeThreshold) {
          ChangeClassification <-
            paste(
              "With phenotype threshold value ",
              phenotypeThreshold,
              " - significant in males (",
              direction_males,
              "), females (",
              direction_females,
              ") and in combined dataset (",
              direction_all,
              ")",
              sep = ""
            )
        }
        # combined & males & !females
        if (all_p.value < phenotypeThreshold &&
          male_p.value < phenotypeThreshold &&
          female_p.value >= phenotypeThreshold) {
          ChangeClassification <-
            paste(
              "With phenotype threshold value ",
              phenotypeThreshold,
              " - significant in males (",
              direction_males,
              ") and in combined dataset (",
              direction_all,
              ")",
              sep = ""
            )
        }
        # combined & !males & females
        if (all_p.value < phenotypeThreshold &&
          male_p.value >= phenotypeThreshold &&
          female_p.value < phenotypeThreshold) {
          ChangeClassification <-
            paste(
              "With phenotype threshold value ",
              phenotypeThreshold,
              " - significant in females (",
              direction_females,
              ") and in combined dataset (",
              direction_all,
              ")",
              sep = ""
            )
        }
        # combined & !males & !females
        if (all_p.value < phenotypeThreshold &&
          male_p.value >= phenotypeThreshold &&
          female_p.value >= phenotypeThreshold) {
          ChangeClassification <-
            paste(
              "With phenotype threshold value ",
              phenotypeThreshold,
              " - significant in combined dataset only (",
              direction_all,
              ")",
              sep = ""
            )
        }
        # !combined & males & females
        if (all_p.value >= phenotypeThreshold &&
          male_p.value < phenotypeThreshold &&
          female_p.value < phenotypeThreshold) {
          ChangeClassification <-
            paste(
              "With phenotype threshold value ",
              phenotypeThreshold,
              " - significant in males (",
              direction_males,
              ") and females (",
              direction_females,
              ") datasets",
              sep = ""
            )
        }
        # !combined & males & !females
        if (all_p.value >= phenotypeThreshold &&
          male_p.value < phenotypeThreshold
        && female_p.value >= phenotypeThreshold) {
          ChangeClassification <-
            paste(
              "With phenotype threshold value ",
              phenotypeThreshold,
              " - significant in males (",
              direction_males,
              ") dataset only",
              sep = ""
            )
        }
        # !combined & !males & females
        if (all_p.value >= phenotypeThreshold &&
          male_p.value >= phenotypeThreshold
        && female_p.value < phenotypeThreshold) {
          ChangeClassification <-
            paste(
              "With phenotype threshold value ",
              phenotypeThreshold,
              " - significant in females (",
              direction_females,
              ") dataset only",
              sep = ""
            )
        }
      } else {
        if (all_p.value < phenotypeThreshold) {
          ChangeClassification <-
            paste0(
              "With phenotype threshold value ",
              phenotypeThreshold,
              " - significant for the sex (",
              lsex,
              ") tested (",
              direction_all,
              ")"
            )
        }
      }
    }
    else {
      ChangeClassification <- NA
    }
  }
  ##############
  if (!is.na(ChangeClassification)) {
    outList <- list(
      "Classification tag"                 = ChangeClassification,
      "Sex in the input model"             = SexInTheCleanedInputModel,
      "Overall p-value threshold"          = phenotypeThreshold,
      "Sex specific p-value threshold"     = SexSpecificThreshold,
      "Active Sex levels"                  = lsex
    )
  } else {
    outList <- NULL
  }
  return(outList)
}
## ------------------------------------------------------------------------------
