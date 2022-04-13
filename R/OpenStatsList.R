setClass(
  "OpenStatsList",
  representation(
    datasetPL = "data.frame",
    refGenotype = "character",
    testGenotype = "character",
    hemiGenotype = "character",
    clean.dataset = "logical",
    dataset.colname.genotype = "character",
    dataset.colname.sex = "character",
    dataset.colname.batch = "character",
    dataset.colname.lifestage = "character",
    dataset.colname.weight = "character",
    dataset.values.missingValue = "character",
    dataset.values.male = "character",
    dataset.values.female = "character",
    dataset.values.early = "character",
    dataset.values.late = "character",
    datasetUNF = "data.frame"
  )
)
OpenStatsList <- function(dataset,
                          testGenotype = "experimental",
                          refGenotype = "control",
                          hemiGenotype = NULL,
                          clean.dataset = TRUE,
                          dataset.colname.genotype = "biological_sample_group",
                          dataset.colname.sex = "sex",
                          dataset.colname.batch = "date_of_experiment",
                          dataset.colname.lifestage = "LifeStage",
                          dataset.colname.weight = "weight",
                          dataset.values.missingValue = c(" ", ""),
                          dataset.values.male = NULL,
                          dataset.values.female = NULL,
                          dataset.values.early = NULL,
                          dataset.values.late = NULL,
                          debug = TRUE) {
  r <- suppressMessagesANDWarnings(OpenStatsList0(
    dataset,
    testGenotype                = testGenotype,
    refGenotype                 = refGenotype,
    hemiGenotype                = hemiGenotype,
    clean.dataset               = clean.dataset,
    dataset.colname.genotype    = dataset.colname.genotype,
    dataset.colname.sex         = dataset.colname.sex,
    dataset.colname.batch       = dataset.colname.batch,
    dataset.colname.lifestage   = dataset.colname.lifestage,
    dataset.colname.weight      = dataset.colname.weight,
    dataset.values.missingValue = dataset.values.missingValue,
    dataset.values.male         = dataset.values.male,
    dataset.values.female       = dataset.values.female,
    dataset.values.early        = dataset.values.early,
    dataset.values.late         = dataset.values.late
  ),
  sup.messages = !debug,
  sup.warnings = FALSE
  )
  return(r)
}


OpenStatsList0 <-
  function(dataset,
           testGenotype = "experimental",
           refGenotype = "control",
           hemiGenotype = NULL,
           clean.dataset = TRUE,
           dataset.colname.genotype = "biological_sample_group",
           dataset.colname.sex = "sex",
           dataset.colname.batch = "date_of_experiment",
           dataset.colname.lifestage = "LifeStage",
           dataset.colname.weight = "weight",
           dataset.values.missingValue = c(" ", ""),
           dataset.values.male = NULL,
           dataset.values.female = NULL,
           dataset.values.early = NULL,
           dataset.values.late = NULL) {
    testGenotype <- as.character(testGenotype)
    refGenotype <- as.character(refGenotype)
    if (is.null(dataset) || !is0(dataset, "data.frame")) {
      message0("error ~> Null dataset or not a data.frame.")
      return(NULL)
    }
    message0(
      "Input data of the dimensions, ",
      paste(
        c("rows", "columns"),
        dim(dataset),
        sep = " = ",
        collapse = ", "
      )
    )
    dataset_unfiltered <- dataset
    if (clean.dataset) {
      sta.time <- Sys.time()
      message0("Checking the input data in progress ...")
      ## Replace missing values specified in the user format with NA
      dataset <- droplevels0(dfNAreplce(df = dataset, NAsymbol = dataset.values.missingValue))
      chkcols <- checkOpenStatsColumns(
        dataset = dataset,
        vars = c(
          dataset.colname.genotype,
          dataset.colname.sex,
          dataset.colname.batch,
          dataset.colname.lifestage,
          dataset.colname.weight
        )
      )
      dataset <- dataset[, order(names(dataset)), drop = FALSE]
      ###############################
      dataset <- fIsBecome(dataset,
        is = dataset.colname.genotype,
        rename = "Genotype",
        ifexist = "Genotype"
      )
      if (!"Genotype" %in% names(dataset) ||
        all(!chkcols) ||
        length(c(testGenotype, refGenotype)) != 2) {
        message0(
          "error ~> Please make sure `dataset.colname.xxx` and/or Genotype column/levels are properly specified."
        )
        return(NULL)
      }
      ###############################
      dataset <- fIsBecome(dataset,
        is = dataset.colname.sex,
        rename = "Sex",
        ifexist = "Sex"
      )
      dataset <- fIsBecome(
        dataset = dataset,
        is = dataset.colname.lifestage,
        rename = "LifeStage",
        ifexist = "LifeStage"
      )
      dataset <- fIsBecome(
        dataset = dataset,
        is = c(dataset.colname.batch),
        rename = "Batch",
        ifexist = "Batch"
      )
      dataset <- fIsBecome(dataset,
        is = dataset.colname.weight,
        rename = "Weight",
        ifexist = "Weight"
      )
      dataset <- factorise(dataset, c("Genotype", "Sex", "Batch", "LifeStage"))
      dataset <- droplevels0(dataset)
      ## Renew levels
      dataset <- OpenStatsListRelabling(
        dataset = dataset,
        col = "Sex",
        l1 = dataset.values.female,
        rel1 = "Female",
        l2 = dataset.values.male,
        rel2 = "Male"
      )
      dataset <- OpenStatsListRelabling(
        dataset = dataset,
        col = "LifeStage",
        l1 = dataset.values.early,
        rel1 = "Early",
        l2 = dataset.values.late,
        rel2 = "Late"
      )
      ## Hemi to test genotype replacement
      if (!is.null(hemiGenotype)) {
        if (any(rownames(dataset[dataset$Genotype == hemiGenotype, ]))) {
          levels(dataset$Genotype)[levels(dataset$Genotype) == hemiGenotype] <-
            testGenotype
          message0(
            "Hemizygotes `",
            hemiGenotype,
            "` have been relabelled to test genotype `",
            testGenotype,
            "`"
          )
        }
      }
      ## Clean genotypes
      if (length(setdiff(
        rownames(dataset),
        rownames(dataset[dataset$Genotype %in% c(testGenotype, refGenotype), ])
      )) >
        0) {
        if (any(!c(testGenotype, refGenotype) %in% levels(dataset$Genotype))) {
          message0(
            "error ~> Mismatch between `Genotype` levels and input levels:",
            "\n\t Genotype Levels:\n\t\t",
            pasteComma(sort0(levels(dataset$Genotype)), sep = ",\n\t\t"),
            "\n\t Input levels   :\n\t\t",
            pasteComma(sort0(c(
              testGenotype, refGenotype
            )), replaceNull = TRUE, sep = "\n\t\t")
          )
          return(NULL)
        }
        dataset <- subset(
          dataset,
          dataset$Genotype %in% c(testGenotype, refGenotype)
        )
        message0(
          "Dataset has been cleaned to only keep the `Genotype` values: \n\t  ",
          pasteComma(numberingX(c(
            testGenotype,
            refGenotype
          )),
          replaceNull = TRUE,
          sep = ",\n\t  "
          )
        )
      }
      # NULL is ignored as it is checked below (avoid double messages)
      if (!is.null(dataset) && is.df.empty(dataset)) {
        message0("Check failed ~> The empty dataset or the preprocessing ended up an empty dataset.")
        return(NULL)
      }
      ## Clean the empty records!
      dataset <- CleanEmptyRecords(dataset, c("Genotype", "Sex", "Batch", "LifeStage"))
      ## CHECKS
      dataset <- checkDataset(
        droplevels(dataset),
        testGenotype,
        refGenotype
      )
      if (is.null(dataset)) {
        return(NULL)
      }

      if ("Weight" %in% colnames(dataset)) {
        if (!is.numeric(dataset$Weight)) {
          message0("`Weight` values are not numeric then renamed to `Weight_labels`")
          colnames(dataset)[colnames(dataset) == "Weight"] <-
            "Weight_labels"
        }
        coNames <- names(dataset) %in% c("Genotype", "Sex", "LifeStage")
        wsglvls <- tapply(X = dataset$Weight, INDEX = interaction(dataset[, coNames]), function(x) {
          length(na.omit(x))
        }, default = 0)
        message0(
          "Total `Weight` data points for ",
          MakeCaptionFromNamesAndIndecis(name = names(dataset), index = coNames),
          ": \n\t Level(frequency): \n\t  ",
          pasteComma(numberingX(paste0(
            names(wsglvls), "(", wsglvls, ")"
          )),
          truncate = FALSE,
          sep = ",\n\t  "
          )
        )
        if (min(wsglvls, na.rm = TRUE) <= 2 &&
          sum(is.na(dataset$Weight)) > 0) {
          message0(
            "`Weight` column has (<2) data points for at least one level of ",
            MakeCaptionFromNamesAndIndecis(name = names(dataset), index = coNames),
            ".\n\t The `Weight` column renamed to `Weight_labels`"
          )
          colnames(dataset)[colnames(dataset) == "Weight"] <-
            "Weight_labels"
        }
      }
      message0(
        "Successfully performed checks in ",
        round(difftime(Sys.time(), sta.time, units = "sec"), 2),
        " second(s)."
      )
    } else {
      message0("No check performed on the input data")
    }
    r <- new(
      "OpenStatsList",
      datasetPL = as.data.frame(dataset),
      refGenotype = as.character(refGenotype),
      testGenotype = as.character(testGenotype),
      hemiGenotype = ifelse(is.null(hemiGenotype), character(0), hemiGenotype),
      dataset.colname.batch = ifelse(
        is.null(dataset.colname.batch),
        character(0),
        dataset.colname.batch
      ),
      dataset.colname.lifestage = ifelse(
        is.null(dataset.colname.lifestage),
        character(0),
        dataset.colname.lifestage
      ),
      dataset.colname.genotype = ifelse(
        is.null(dataset.colname.genotype),
        character(0),
        dataset.colname.genotype
      ),
      dataset.colname.sex = ifelse(
        is.null(dataset.colname.sex),
        character(0),
        dataset.colname.sex
      ),
      dataset.colname.weight = ifelse(
        is.null(dataset.colname.weight),
        character(0),
        dataset.colname.weight
      ),
      dataset.values.missingValue = ifelse(
        is.null(dataset.values.missingValue),
        character(0),
        dataset.values.missingValue
      ),
      dataset.values.male = ifelse(
        is.null(dataset.values.male),
        character(0),
        as.character(dataset.values.male)
      ),
      dataset.values.female = ifelse(
        is.null(dataset.values.female),
        character(0),
        as.character(dataset.values.female)
      ),
      dataset.values.early = ifelse(
        is.null(dataset.values.early),
        character(0),
        as.character(dataset.values.early)
      ),
      dataset.values.late = ifelse(
        is.null(dataset.values.late),
        character(0),
        as.character(dataset.values.late)
      ),
      clean.dataset = as.logical(clean.dataset),
      datasetUNF = as.data.frame(dataset_unfiltered)
    )
    return(invisible(r))
  }


#-------------------------------------------------------------------------------
## Check dataset for the minimum required info and additional cleaning steps
checkDataset <- function(dataset,
                         testGenotype,
                         refGenotype = "+/+") {
  dataset <- droplevels(dataset)
  if (any(c("Genotype", "Sex") %in% colnames(dataset))) {
    coNames <- names(dataset) %in% c("Genotype", "Sex", "LifeStage")
    InGS <- interaction(dataset[, coNames])
    tbGSL <- table(InGS)


    message0(
      "Total samples in ",
      MakeCaptionFromNamesAndIndecis(name = names(dataset), index = coNames),
      ": \n\t Level(frequency): \n\t  ",
      pasteComma(
        numberingX(paste0(names(tbGSL), "(", tbGSL, ")")),
        truncate = FALSE,
        sep = "\n\t  "
      )
    )
    if (min(tbGSL) < 1) {
      message0(
        "No observations detected in ",
        MakeCaptionFromNamesAndIndecis(name = names(dataset), index = coNames),
        " for:\n\t",
        pasteComma(names(tbGSL[tbGSL < 1]), truncate = FALSE, sep = ",\n\t")
      )
    }
    dataset <- droplevels(dataset[InGS %in% names(tbGSL[tbGSL >= 1]), , drop = FALSE])
    ## Check of genotype and sex levels after cleaning
    if (nlevels(dataset$Genotype) != 2) {
      message0(
        "error ~> `Genotype` column must have two levels. Current levels:\n\t ",
        pasteComma(levels(dataset$Genotype), sep = ",\n\t")
      )
      return(NULL)
    }
    if ("Sex" %in% names(dataset) &&
      nlevels(dataset$Sex) > 2) {
      message0(
        "error ~> `Sex` column must have one or two levels. Current levels:\n\t ",
        pasteComma(levels(dataset$Sex), sep = ",\n\t")
      )
      return(NULL)
    }
    ## Check for sex levels - we want to have 'Female' and/or 'Male' only
    wrong_sex_levels <- setdiff(levels(dataset$Sex), c("Female", "Male"))
    if ("Sex" %in% names(dataset) &&
      !length(wrong_sex_levels) < 1) {
      message0(
        "error ~> Sex has undefined levels. See:\n\t ",
        pasteComma(wrong_sex_levels, truncate = FALSE, sep = ",\n\t")
      )
      return(NULL)
    }
    ## Check for reference genotype records
    if (refGenotype %in% levels(dataset$Genotype)) {
      dataset$Genotype <- relevel(dataset$Genotype, ref = refGenotype)
    }
  } else {
    message0(
      "error ~ Neither `Genotype` or `Sex` found in the input data."
    )
    return(NULL)
  }
  return(dataset)
}

OpenStatsListBuilder <- function(PhenListobject,
                                 DOE = NULL,
                                 DOB = NULL,
                                 d.threshold = 16 * 7,
                                 debug = TRUE) {
  #### Negative age will be removed
  PhenListobject@datasetPL <- droplevels(PhenListobject@datasetPL)

  Ageing <- !is.null(DOE) &&
    !is.null(DOB) &&
    all(c(DOE, DOB) %in% names(PhenListobject@datasetPL))

  if (Ageing) {
    age.in.day <- as.Date(PhenListobject@datasetPL[, DOE]) - as.Date(PhenListobject@datasetPL[, DOB])
    PhenListobject@datasetPL$Age <- as.numeric(age.in.day)
    PhenListobject@datasetPL$LifeStage <- ifelse(age.in.day > d.threshold, "Late", "Early")
    PhenListobject@datasetPL$LifeStage <- as.factor(PhenListobject@datasetPL$LifeStage)
    PhenListobject@datasetPL <- PhenListobject@datasetPL[PhenListobject@datasetPL$Age > 0, ]
    message0("Age range: ", paste0(range(age.in.day), collapse = "-"))
  } else {
    message0("DOE and DOB are not specified. Then PhenList is returned.")
  }

  LL <- levels(PhenListobject@datasetPL$LifeStage)
  LS <- levels(PhenListobject@datasetPL$Sex)
  LG <- levels(PhenListobject@datasetPL$Genotype)

  if (length(LL) != 2 && debug && Ageing) {
    message0(
      "Ageing pipeline requires two levels in the LifeStage. Levels: ",
      pasteComma(LL,
        replaceNull = FALSE,
        truncate = FALSE,
        sep = ",\n\t\t"
      ),
      "\nNormal pipeline will apply to this data"
    )
    # return(NULL)
  }
  if (length(LS) != 2 && debug) {
    message0(
      "There should be two levels in Sex. Levels:\n\t ",
      pasteComma(LS,
        replaceNull = FALSE,
        truncate = FALSE,
        sep = ",\n\t"
      )
    )
    # return(NULL)
  }
  if ("Weight" %in% names(PhenListobject@datasetPL) &&
    sum(is.na(PhenListobject@datasetPL[, "Weight"])) > 0 && debug) {
    message0("There are ", sum(is.na(PhenListobject@datasetPL[, "Weight"])), " NAs in body weights.")
  }

  if (length(LG) < 2 && debug) {
    message0(
      "Check failed ~> Genotype must have two levels. Levels: \n\t",
      pasteComma(LG,
        replaceNull = FALSE,
        truncate = FALSE,
        sep = ",\n\t"
      )
    )
    return(NULL)
  }
  OpenStatsList <- unclass(PhenListobject)
  class(OpenStatsList) <- "OpenStatsList"
  return(OpenStatsList)
}

summary.OpenStatsList <- function(object,
                                  vars = NULL,
                                  ...) {
  requireNamespace("summarytools")
  r <- NULL
  df <- SelectVariablesOrDefault(data = object@datasetPL, vars)
  if (!is.null(df)) {
    message0("Working on the summary table ...")
    r <- summarytools::dfSummary(df, justify = "l", ...)
  }
  return(r)
}

print.OpenStatsList <- function(x,
                                vars = NULL,
                                ...) {
  r <- summary.OpenStatsList(object = x, vars = vars, ...)
  return(invisible(r))
}

plot.OpenStatsList <- function(x,
                               vars = NULL,
                               ...) {
  requireNamespace("Hmisc")
  df <- SelectVariablesOrDefault(data = x@datasetPL, vars)
  if (!is.null(df)) {
    message0("Working on the plot ...")
    r <- Hmisc::describe(df)
    plot <- suppressWarnings(plot(r, ...))
    return(plot)
  } else {
    message0("No variable found in the data. Please make sure that `vars` exist in the input data")
    return(NULL)
  }
}

SelectVariablesOrDefault <- function(data, vars = NULL) {
  cnames <- c(
    "Genotype",
    "Sex",
    "LifeStage",
    "Batch",
    "age_in_weeks",
    "phenotyping_center",
    "metadata_group"
  )
  ###################
  cnamesF <- varExistsInDF(data = data, if (is.null(vars)) {
    cnames
  } else {
    vars
  })
  if (is.null(cnamesF) || length(cnamesF) < 1) {
    message0(
      "Non of the specified variables exist in the data. See the list below:\n\t",
      pasteComma(cnames, truncate = FALSE, sep = ",\n\t")
    )
    return(NULL)
  }
  return(data[, cnamesF, drop = FALSE])
}

varExistsInDF <- function(data = NULL, vars = NULL) {
  if (is.null(vars) || is.null(data)) {
    return(NULL)
  }
  r <- vars[vars %in% names(data)]
  return(r)
}

MakeCaptionFromNamesAndIndecis <- function(name, index) {
  r <- if (length(name[index]) > 1) {
    paste0(paste(name[index], collapse = ":", sep = ":"), " interactions")
  } else {
    name[index]
  }
  return(r)
}
