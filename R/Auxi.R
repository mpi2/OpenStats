# Startup message
.onAttach <- function(lib, pkg) {
  packageStartupMessage(
    paste0(
      "\n >===============================================================================<",
      "\n OpenStats is developed by International Mouse Phenotyping Consortium (IMPC) ",
      "\n More details           : https://www.mousephenotype.org/                         ",
      "\n Source code and issues : https://git.io/Jv5w0                                    ",
      "\n Contact us             : hamedhm@ebi.ac.uk                                       ",
      "\n >===============================================================================<"
    ),
    domain = NULL,
    appendLF = TRUE
  )
}

asFactorAndSelectVariable <- function(x = NULL, col = NULL) {
  x <- droplevels0(x)
  if (!is.null(x) &&
    !is.null(col) &&
    col %in% names(x)) {
    if (length(col) > 1) {
      message0("Only one element of the `col` parameter will be used.")
    }
    r <- as.factor(x[, col[1]])
  } else {
    r <- NULL
  }
  return(r)
}

ML2REML <- function(x, debug = FALSE) {
  if (!is.null(x) &&
    is(x, c("lme", "gls"))) {
    if (debug) {
      message0("\tRecovering the REML object from ML ...")
    }
    x <- tryCatch(
      expr = update(x, method = "REML"),
      error = function(e) {
        message0("\t\tError(s) in updating ML object to REML. See: ")
        message0("\t\t", e, breakLine = FALSE)
        return(NULL)
      },
      warning = function(w) {
        message0("\t\tWarning(s) in updating ML object to REML. See: ")
        message0("\t\t", w, breakLine = FALSE)
        return(NULL)
      }
    )
  }
  return(x)
}

REML2ML <- function(x, debug = FALSE) {
  if (!is.null(x) &&
    is(x, c("lme", "gls"))) {
    if (debug) {
      message0("\tCoverting REML object to ML ...")
    }
    x <- tryCatch(
      expr = update(x, method = "ML"),
      error = function(e) {
        message0("\t\tError(s) in updating REML object to ML. See: ")
        message0("\t\t", e, breakLine = FALSE)
        return(NULL)
      },
      warning = function(w) {
        message0("\t\tWarning(s) in updating REML object to ML. See: ")
        message0("\t\t", w, breakLine = FALSE)
        return(NULL)
      }
    )
  }
  return(x)
}



Matrix2List <- function(x, ...) {
  if (is.null(x)) {
    return(NULL)
  }

  if (length(x) == 1 || is(x, "numeric")) {
    return(as.list(x))
  }

  if (!is(x, "matrix")) {
    x <- as.matrix(x)
  }
  r <- as.list(unmatrix0(x, ...))
  return(r)
}

replaceNull <- function(x, replaceBy = "NULL") {
  x <- sapply(x, function(xx) {
    if (is.null(xx)) {
      replaceBy
    } else {
      xx
    }
  })
  return(unlist(x))
}
pasteComma <- function(...,
                       replaceNull = TRUE,
                       truncate = TRUE,
                       width = 100,
                       trailingSpace = TRUE,
                       replaceNullby = "NULL",
                       sep = ",") {
  sep <- ifelse(trailingSpace && sep %in% ",", paste0(sep, " "), sep)
  if (replaceNull) {
    r <- paste(
      replaceNull(list(...), replaceBy = replaceNullby),
      sep      = sep,
      collapse = sep
    )
  } else {
    r <- paste(...,
      sep      = sep,
      collapse = sep
    )
  }

  if (truncate) {
    r <- truncate_text(r, width)
  }

  return(r)
}

truncate_text <- function(x, width) {
  ifelse(nchar(x) > width, paste0(strtrim(x, width), " ..."), x)
}

pasteUnderscore <- function(...) {
  paste(..., sep = "_", collapse = "_")
}
pastedot <- function(...) {
  paste(..., sep = ".", collapse = ".")
}
checkModelTermsInData <- function(formula,
                                  data,
                                  responseIsTheFirst = TRUE,
                                  pattern = "[.~+-()]") {
  formula <- as.formula(formula)
  vars <- all_vars0(formula, functions = FALSE)
  vars <- vars[!grepl(
    pattern = pattern,
    x = vars,
    fixed = TRUE
  )]
  if (responseIsTheFirst) {
    if (!(vars[1] %in% names(data))) {
      message0(
        "Response does not exist in the data!\n\tFormula: ",
        printformula(formula)
      )
      stop(
        "Response has not been specified properly. Please check that the response exists in the data"
      )
    }
  }
  In <- vars %in% names(data)
  if (any(!In)) {
    message0(
      "Some terms in the model are not included in the data. See: \n\t  ",
      pasteComma(vars[!In], replaceNull = FALSE, truncate = FALSE),
      "\n\t Initial  model: ",
      printformula(formula)
    )
    ft <- vars [!In]
    formula <- update.formula(
      formula,
      reformulate0(
        termlabels = c(".", ft, paste0(ft, ":.")),
        response = NULL,
        intercept = TRUE,
        sep = "-"
      )
    )
    message0("\t Polished model: ", printformula(formula))
  }
  return(formula)
}




dfNAreplce <- function(df,
                       NAsymbol = NULL,
                       replceNaBy = NA) {
  if (!is.null(NAsymbol) &&
    length(NAsymbol) > 0) {
    message0(
      "Checking the specified missing values [x",
      length(NAsymbol),
      "] (",
      pasteComma(paste0("`", NAsymbol, "`")),
      ") ..."
    )
    #### This is much faster than df %in% NAsymbol
    n <- length(NAsymbol)
    for (i in seq_along(NAsymbol)) {
      NAs <- NAsymbol[i]
      message0("\t", i, "/", n, ". Checking (`", NAs, "`) ...")
      df[df == NAs] <- replceNaBy
    }
  }
  return(df)
}

reformulate0 <- function(termlabels,
                         response = NULL,
                         intercept = TRUE,
                         sep = "+") {
  if (!is.character(termlabels) || !length(termlabels)) {
    stop("'termlabels' must be a character vector of length at least one")
  }
  has.resp <- !is.null(response)
  termtext <- paste(if (has.resp) {
    "response"
  },
  "~",
  paste(termlabels, collapse = sep),
  collapse = ""
  )
  if (!intercept) {
    termtext <- paste(termtext, "- 1")
  }
  rval <- eval(parse(text = termtext, keep.source = FALSE)[[1L]])
  if (has.resp) {
    rval[[2L]] <- if (is.character(response)) {
      as.symbol(response)
    } else {
      response
    }
  }
  environment(rval) <- parent.frame()
  rval
}

suppressMessagesANDWarnings <- function(exp,
                                        sup.messages = TRUE,
                                        sup.warnings = FALSE) {
  if (sup.messages && sup.warnings) {
    suppressMessages(suppressWarnings(exp))
  } else if (sup.messages && !sup.warnings) {
    suppressMessages(exp)
  } else if (!sup.messages && sup.warnings) {
    suppressWarnings(exp)
  } else {
    exp
  }
}
# Typical fixed effect
TypicalModel <- function(depVariable,
                         withWeight = TRUE,
                         Sex = TRUE,
                         LifeStage = TRUE,
                         mandatory = "Genotype",
                         data,
                         others = NULL,
                         debug = TRUE) {
  colNames <- colnames(data)
  if (!mandatory %in% colNames) {
    stop("Genotype does not found in the dataset!")
  }

  fixed <- reformulate(termlabels = mandatory, response = depVariable)
  if (Sex && ("Sex" %in% colNames)) {
    fixed <- update(fixed, ~ . * Sex)
  }
  if (LifeStage && ("LifeStage" %in% colNames)) {
    fixed <- update(fixed, ~ . * LifeStage)
  }
  if (withWeight && ("Weight" %in% colNames)) {
    fixed <- update(fixed, ~ . + Weight)
  }
  if (!is.null(others)) {
    fixed <- update(fixed, reformulate(response = NULL, termlabels = c(".", others)))
  }
  ####
  if (debug) {
    message0("Initial model: ", printformula(fixed))
  }
  return(fixed)
}



FeasibleTermsInContFormula <- function(formula, data) {
  if (is.null(formula) || is.null(data)) {
    message0("Null data or the formula. Check the data and/or formula")
    stop()
  }
  Allvars <- all_vars0(formula)[all_vars0(formula) %in% names(data)]
  isCat <- !is.continuous(data[, Allvars, drop = FALSE])
  vars <- Allvars[isCat]
  lvars <- length(vars)
  names <- r <- NULL
  if (getResponseFromFormula(formula = formula) %in% vars) {
    message0("\t Response is included in the checks ....")
  }
  if (lvars > 0) {
    for (i in seq_len(lvars)) {
      message0(
        "\t",
        i,
        " of ",
        lvars,
        ". Checking for the feasibility of terms and interactions ..."
      )
      cmb <- combn(vars, i)
      for (j in seq_len(ncol(cmb))) {
        message0("\t\t Checking ", pasteComma(cmb[, j]), " ...")
        xtb <- xtabs(
          formula = paste0("~", paste0(cmb[, j], collapse = "+")),
          data = data,
          drop.unused.levels = FALSE
        )
        r <- c(r, if (all(dim(xtb) >= 2)) {
          min(xtb, na.rm = TRUE)
        } else {
          0
        })
        names <- c(names, paste0(cmb[, j], collapse = ":"))
      }
    }
    return(data.frame(
      names = names,
      min.freq = r,
      stringsAsFactors = FALSE
    ))
  } else {
    return(NULL)
  }
}

variablesInData <- function(df, names, debug = TRUE) {
  if (is.null(df) || is.null(names) || sum(names %in% names(df)) < 1) {
    return(NULL)
  }
  newNames <- names[names %in% names(df)]
  if (debug) {
    message0(
      "Variables that being found in data: ",
      pasteComma(newNames, truncate = FALSE)
    )
  }
  return(newNames)
}

ComplementaryFeasibleTermsInContFormula <- function(formula, data) {
  message0(
    "Checking for the feasibility of terms and interactions ...\n\t Formula: ",
    printformula(formula)
  )
  fbm <- FeasibleTermsInContFormula(formula = formula, data = data)
  if (!is.null(fbm) &&
    (min(fbm$min.freq, na.rm = TRUE) < 1 ||
      length(formulaTerms(formula)) != nrow(fbm))) {
    formula <- update.formula(
      old = formula,
      new =
        reformulate0(
          termlabels = c(".", fbm$names[fbm$min.freq <= 0]),
          response = ".",
          intercept = TRUE,
          sep = "-"
        )
    )
    if (min(fbm$min.freq, na.rm = TRUE) < 1) {
      message0(
        'The following term(s) removed because there is either "no data" or "no data for the interactions":\n\t ** Note. Not all terms necessarily in the initial model \n\t ',
        pasteComma(fbm[fbm$min.freq <= 0, c("names")], replaceNull = FALSE, truncate = FALSE)
      )
    }
  }
  return(formula)
}

sign0 <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  if (sign(x) > 0) {
    return("positive")
  } else if (sign(x) == 0) {
    return("neutral")
  } else if (sign(x) < 0) {
    return("negative")
  } else {
    return(NULL)
  }
}

dist0 <- function(x, func = lower.tri) {
  if (is.null(x)) {
    return(NULL)
  }
  out <- outer(x, x, `-`)
  r <- out[func(out)]
  return(r)
}

CheckMissing <- function(data, formula) {
  if (is.null(formula) || is.null(data)) {
    message0("Null data or the formula. Check the data and/or formula")
    stop()
  }
  org.data <- data
  new.data <- data[complete.cases(data[, all_vars0(formula)]), ]
  missings <- ifelse(all(dim(org.data) == dim(new.data)), 0, dim(org.data)[1] -
    dim(new.data)[1])
  if (missings) {
    message0(
      "The data (variable(s) = ",
      pasteComma(all_vars0(formula), truncate = FALSE),
      ") contain ",
      missings,
      " missing(s) ...\n\tMissing data removed"
    )
  }
  return(invisible(
    list(
      org.data = org.data,
      new.data = droplevels0(new.data),
      missings = missings
    )
  ))
}

droplevels0 <- function(x, ...) {
  if (is.null(x) ||
    class(x) %in% c("matrix", "integer", "double", "numeric")) {
    return(x)
  }
  return(droplevels(x, ...))
}

range0 <- function(x, ...) {
  ran <- range(x, na.rm = TRUE)
  if (length(ran) == 2) {
    return(diff(ran))
  } else {
    return(NULL)
  }
}

order0 <- function(x, levels = FALSE) {
  if (is.null(x)) {
    return(NULL)
  }

  if (levels) {
    r <- x[order(levels(x))]
  } else {
    r <- x[order(x)]
  }
  return(r)
}


FormulaContainsFunction <- function(formula) {
  message0("Checking the input model for including functions ...")
  if (is.null(formula)) {
    message0("Blank formula!")
    return(NULL)
  }
  #
  fFull <- all_vars0(formula, functions = TRUE)
  FAbstract <- all_vars0(formula, functions = FALSE)
  r <- !identical(
    fFull    [grepl(pattern = "[0-9A-Za-z^\\/]", x = fFull)],
    FAbstract[grepl(pattern = "[0-9A-Za-z]", x = FAbstract)]
  )
  if (r) {
    message0("\t Function detected in the input model ...")
  }

  return(r)
}

percentageChangeCont <- function(model,
                                 data,
                                 variable,
                                 depVar,
                                 individual = TRUE,
                                 mainEffsOnlyWhenIndivi = "Sex",
                                 FUN = range0,
                                 sep = " ") {
  if (!(
    !is.null(data) &&
      (!is.null(variable) || !is.null(mainEffsOnlyWhenIndivi)) &&
      !is.null(model) &&
      !is.null(FUN(data[, depVar])) &&
      FUN(data[, depVar]) != 0
  )) {
    return(NULL)
  }
  ####
  model <-
    tryCatch(
      expr = update(model,
        data = NormaliseDataFrame(data)
      ),
      error = function(e) {
        message0(
          "\t\tError(s) in the (combined) effect size estimation for",
          pasteComma(variable),
          ". See: "
        )
        message0("\t\t", e, breakLine = FALSE)
        return(NULL)
      },
      warning = function(w) {
        message0(
          "\t\tWarning(s) in the (combined) effect size estimation for",
          pasteComma(variable),
          ". See: "
        )
        message0("\t\t", w, breakLine = FALSE)
        return(NULL)
      }
    )
  if (is.null(model)) {
    return(NULL)
  }
  coefs <- unlist(
    suppressMessagesANDWarnings(
      summary(model, verbose = FALSE)$tTable[, 1],
      sup.messages = TRUE,
      sup.warnings = TRUE
    )
  )
  ######
  ran <- 1 # --ONE-- Just to cancel 100 but keep the percent (*100) FUN(data[, depVar])
  if (is.null(coefs) || is.null(ran) || is.nan(ran)) {
    return(NULL)
  }
  ######
  if (individual) {
    out <- sapply(variable, function(x) {
      message0(
        "\tCalculating the percentage change for: ",
        pasteComma(x)
      )
      if (is.numeric(data[, x])) {
        r <- coefs[-1]
        names(r) <- NULL
      } else {
        r <- coefs
        lr <- length(r)
        ol <- order0(levels(data[, x]))
        ##########
        if (lr != nlevels(data[, x])) {
          message0(
            "\tCare reguired for the percentage change calculation for: ",
            x,
            ". Levels = ",
            names(r),
            "Coefficients: ",
            r
          )
          names(r) <- paste(ol[seq_len(lr)], "_CareRequired")
        } else {
          names(r) <- ol
        }
      }
      return(r / ran * 100)
    }, USE.NAMES = TRUE)
    return(Matrix2List(out, sep = sep))
  } else {
    out <- extractCoefOfInterest(
      coefs = coefs,
      main = mainEffsOnlyWhenIndivi,
      data = data
    )
    return(
      list(
        "Value"               = as.list(out),
        "Variable"            = mainEffsOnlyWhenIndivi,
        "Model"               = printformula(formula(model)),
        "Type"                = "Standardized interaction coefficients ",
        "Percentage change"   = Matrix2List(out / ran * 100, sep = sep)
      )
    )
  }
}

optimM <- function(optimise) {
  paste0(c("Fixed term = ", "Weight term = ", "Random effect term = "),
    optimise,
    collapse = ", "
  )
}

applyFormulaToData <- function(formula = NULL, data, add = FALSE) {
  if (is.null(formula) || is.null(all_vars0(formula))) {
    return(data)
  }
  if (is.null(data)) {
    return(NULL)
  }
  nms <-
    trimws(scan(
      text = paste(unlist(as.list(
        attr(terms(as.formula(formula)), "variables")
      ))[-1], sep = "+", collapse = " + "),
      what = "",
      sep = "+",
      quiet = TRUE
    ))

  m <- sapply(nms, function(x) {
    eval(parse(text = x), data)
  })
  if (!is.null(m) && add) {
    m <- cbind(data, m)
  }
  return(list(data = m, names = nms))
}

optimiseMessage <- function(optimise) {
  message0(
    "The model optimisation is ",
    ifelse(
      all(optimise),
      "in progress ...",
      paste0("in the following order:\n\t", optimM(optimise))
    )
  )
}

extractCoefOfInterest <- function(coefs, main = "Sex", data) {
  lvls <- lapply(
    main,
    FUN = function(x) {
      if (!is.numeric(data[, x])) {
        r <- paste(x, levels(data[, x]), sep = "")
        r <- c(paste0(r, ":"), paste0(":", r))
      } else {
        r <- x
      }
      return(r)
    }
  )

  Cnames <- names(coefs)
  r <- coefs[grepl(
    pattern = pasteBracket(
      unlist(lvls),
      left = "",
      right = "",
      col = "|"
    ),
    x = Cnames
  )]
  return(r)
}

pasteBracket <- function(...,
                         replaceNull = TRUE,
                         col = "|",
                         right = "]:",
                         left = "[") {
  if (replaceNull) {
    paste(
      left,
      replaceNull(list(...), replaceBy = "NULL"),
      right,
      sep = "",
      collapse = col
    )
  } else {
    paste(left, ..., right, sep = " ", collapse = col)
  }
}

eff.size <- function(object,
                     data = NULL,
                     depVariable = "data_point",
                     effOfInd = "Genotype",
                     errorReturn = NULL,
                     debug = FALSE) {
  if (all(is.null(data))) {
    data <- NormaliseDataFrame(getData(object))
  }
  f <- reformulate(termlabels = effOfInd, depVariable)
  # Remove NAs
  data <- data[complete.cases(data[, all_vars0(f)]), ]
  if (!(depVariable %in% names(data) &&
    length(na.omit(data[, depVariable])) > 1)) {
    return(NULL)
  }

  agr <- aggregate(f, data = data, FUN = function(x) {
    mean(x, na.rm = TRUE)
  })
  if (debug) {
    cat("\n\n")
    print(agr)
    print(dim(agr))
  }
  if (any(dim(agr) < 2)) {
    message0(
      " \t\t(standardized) Effect size estimation: No variation or less than two levels in ",
      pasteComma(effOfInd, collapse = ",")
    )
    return(errorReturn)
  }

  NModel <-
    tryCatch(
      expr = update(
        object,
        reformulate(
          termlabels = effOfInd,
          response = ".",
          intercept = TRUE
        ),
        data = data
      ),
      error = function(e) {
        message0(
          "\t\tError(s) in the effect size estimation for",
          pasteComma(effOfInd),
          ". See: "
        )
        message0("\t\t", e, breakLine = FALSE)
        return(NULL)
      },
      warning = function(w) {
        message0(
          "\t\tWarning(s) in the effect size estimation for",
          pasteComma(effOfInd),
          ". See: "
        )
        message0("\t\t", w, breakLine = FALSE)
        return(NULL)
      }
    )
  if (!is.null(NModel)) {
    CoefEffSizes <- is.continuous(data[, effOfInd, drop = FALSE])
    PerChange <- percentageChangeCont(
      model = NModel,
      data = data,
      variable = effOfInd,
      depVar = depVariable
    )
    if (sum(CoefEffSizes)) {
      # For continues covariates it is the coefficient
      efSi <- list(
        "Value"               = as.list(coef(NModel))[[effOfInd]][1],
        "Variable"            = effOfInd,
        "Model"               = printformula(formula(NModel)),
        "Type"                = "Standardized coefficient",
        "Percentage change"   = PerChange
      )
    } else {
      # For categorical covariates it is the mean difference
      MDiff <- max(dist(agr[, depVariable, drop = FALSE],
        method = "maximum"
      ),
      na.rm = TRUE
      )
      r <- resid(NModel)
      sd <- sd0(r, na.rm = TRUE)
      efSi <- list(
        "Value" = ifelse(!is.na(sd) &&
          sd > 0, abs(MDiff) / sd, NA),
        "Variable" = effOfInd,
        "Model" = printformula(formula(NModel)),
        "Type" = "Mean differences",
        "Percentage change" = PerChange
      )
    }
  } else {
    efSi <- errorReturn
  }
  return(efSi)
}

noVariation <- function(data, f = "~Genotype") {
  xtb <- xtabs(
    formula = f,
    drop.unused.levels = TRUE,
    data = data
  )
  if (any(dim(xtb) < 2)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

ConvDf2Flat <- function(dframe,
                        ch1 = "*",
                        ch2 = ":",
                        chend = ";") {
  out <- apply(as.data.frame(dframe), 1, function(x) {
    if (length(x) > 2) {
      paste(paste(paste(x[seq_len(length(x) - 1)], collapse = ch1),
        trimws(x[length(x)]),
        sep = ch2
      ),
      collapse = chend
      )
    } else if (length(x) == 2) {
      paste(paste(x[1], ch2, x[2]), collapse = chend)
    } else {
      paste(paste(x), collapse = chend)
    }
  })
  return(out)
}

#### List of all procedures
lop <- function() {
  return(
    c(
      "TYPICAL",
      "ABR",
      "ACS",
      "ALZ",
      "BLK",
      "BWT",
      "CAL",
      "CBC",
      "CHL",
      "CSD",
      "DXA",
      "ECG",
      "ECH",
      "ELZ",
      "EVL",
      "EVM",
      "EVO",
      "EVP",
      "EYE",
      "FER",
      "GEL",
      "GEM",
      "GEO",
      "GEP",
      "GPL",
      "GPM",
      "GPO",
      "GRS",
      "HEM",
      "HIS",
      "HWT",
      "IMM",
      "INS",
      "IPG",
      "OFD",
      "PAT",
      "VIA",
      "XRY"
    )
  )
}

diff0 <- function(x) {
  r <- if (length(x) < 2) {
    x
  } else {
    diff(x)
  }
  return(r)
}

summary1 <- function(x, ...) {
  r <- tryCatch(
    expr = summary(x, ...),
    warning = function(war) {
      message0("\t * The summary failed with the warning (see below): ")
      message0("\t", war, breakLine = FALSE)
      return(NULL)
    },
    error = function(err) {
      message0("\t * The summary failed with the error (see below): ")
      message0("\t   ", err, breakLine = FALSE)
      return(NULL)
    }
  )
  return(r)
}

anova0 <- function(x, ...) {
  r <- tryCatch(
    expr = anova(x, ...),
    warning = function(war) {
      message0("\t * The ANOVA failed with the warning (see below): ")
      message0("\t", war, breakLine = FALSE)
      return(NULL)
    },
    error = function(err) {
      message0("\t * The ANOVA failed with the error (see below): ")
      message0("\t   ", err, breakLine = FALSE)
      return(NULL)
    }
  )
  return(r)
}

summary0 <- function(x, ...) {
  if (is.null(x) || length(x) < 1) {
    message0("Null column found in the data!")
    return(x)
  }
  if (is.numeric(x)) {
    r <- summary1(diff0(x), ...)
  } else {
    xx <- as.factor(x)
    r <- c(sort(levels(xx)), summary1(diff0(as.integer(xx)), ...))
  }
  return(r)
}

RemoveDuplicatedColumnsFromDfandTrimWhiteSpace <- function(x, formula = NULL, trimWS = TRUE) {
  x <- as.data.frame(x)
  if (!is.null(formula) ||
    sum(all_vars0(formula) %in% names(x)) > 1) {
    vars <- all_vars0(formula)
  } else {
    vars <- names(x)
  }
  colVars <- names(x) %in% vars
  if (sum(colVars)) {
    subX <- x[, colVars, drop = FALSE]
    #####################
    if (trimWS) {
      subX <- trimColsInDf(df = subX)
    }
    #####################
    message0(
      "Checking duplications in the data model:\n\t ",
      pasteComma(names(x)[colVars],
        truncate = FALSE
      )
    )

    # numCols = is.continuous(subX)
    # ConCols = subX[,  numCols, drop = FALSE]
    # CatCols = subX[, !numCols, drop = FALSE]
    dcols <- duplicated(lapply(subX, summary0))
    if (any(dcols)) {
      message0(
        "\tDuplicated columns found (and removed) in the input data. Removed variables:\n\t ",
        pasteComma(names(subX)[dcols], truncate = FALSE)
      )
    } else {
      message0("\tNo duplicate found.")
    }
    uniqCols <- subX   [, !dcols, drop = FALSE]
    r <- cbind(x[, !colVars, drop = FALSE], uniqCols)
    return(r)
  } else {
    message0("Formula terms do not exist in the input data")
    return(x)
  }
}

colExists <- function(name, data) {
  if ((name %in% names(data)) &&
    length(complete.cases(data[, name])) > 0) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

listFun <- function(list,
                    FUN,
                    debug = FALSE,
                    messageUnusedArgs = FALSE) {
  if (debug) {
    message0("\tApplied model: ", FUN)
  }
  fArgs <- names(list) %in% formalArgs(args(FUN))
  if (messageUnusedArgs &&
    debug &&
    !identical(list[names(list)[fArgs]], list)) {
    message0(
      "Unused argument(s) in the function input. See:\n\t",
      pasteComma(names(list)[!fArgs],
        truncate = TRUE,
        width = 75
      )
    )
  }
  l <- list[names(list)[fArgs]]
  return(l)
}

RandomEffectCheck <- function(formula, data) {
  if (is.null(formula) || is.null(data)) {
    return(NULL)
  }
  message0(
    "Checking the random effect term ...\n\tFormula: ",
    printformula(formula)
  )
  difTerms <- setdiff(all_vars0(formula), names(data))
  if (length(difTerms)) {
    message0(
      "\tSome terms in the random effect do not exist in the data. See:\n\t ",
      difTerms,
      "\n\tRandom effect is set to NULL"
    )
    return(NULL)
  } else {
    return(formula)
  }
}

ModelChecks <- function(fixed,
                        data,
                        checks = c(0, 0, 0, 0),
                        responseIsTheFirst = TRUE) {
  if (length(checks) != 4) {
    message0('"checks" must be a vector of 4 0/1 TRUE/FALSE values. Example c(1,1,1,1) or c(1,1,1,0) or c(0,0,0,0)')
    return(fixed)
  }
  if (any(checks > 0)) {
    if (checks[1]) {
      fixed <- checkModelTermsInData(
        formula = fixed,
        data = data,
        responseIsTheFirst = TRUE
      )
    }
    if (checks[2]) {
      fixed <- removeSingleLevelFactors(formula = fixed, data = data)
    }
    if (checks[3]) {
      fixed <- removeSingleValueContinuousVariables(formula = fixed, data = data)
    }
    if (checks[4]) {
      fixed <- ComplementaryFeasibleTermsInContFormula(formula = fixed, data = data)
    }
    message0("\tChecked model: ", printformula(fixed))
  }
  fixed <- missingInVariable(
    fixed = fixed,
    data = data,
    threshold = 50
  )
  return(fixed)
}

missingInVariable <- function(fixed = NULL,
                              data = NULL,
                              threshold = 50) {
  message0("Check missings in progress ...")
  if (is.null(data) || nrow(data) < 1) {
    message0("\tNo data at all")
    return(fixed)
  }
  if (is.null(fixed) || length(all_vars0(fixed)) < 1) {
    message0("\tNo formula imported")
    return(fixed)
  }
  vars <- all_vars0(fixed)
  if (length(vars) < 1 || all(!vars %in% names(data))) {
    message0("\tVariables in the formula do not exist in the input data")
    return(fixed)
  }
  newVars <- vars[vars %in% names(data)]
  lapply(newVars, function(v) {
    MissingCounterFunction(
      v = v,
      data = data,
      threshold = threshold
    )
  })
  return(fixed)
}

MissingCounterFunction <- function(v, data, threshold) {
  if (length(data[, v]) < 1) {
    return(NULL)
  }
  missingPercentage <- round(sum(is.na(data[, v])) / length(data[, v]) * 100, 2)
  message0(
    "\t Missings in variable `",
    v,
    "`: ",
    missingPercentage,
    ifelse(
      length(threshold) > 0 &&
        missingPercentage > threshold,
      paste0("% [more than ", threshold, "% of the data]"),
      "%"
    )
  )
}

termInTheModel <- function(model, term, message = FALSE) {
  if (message) {
    message0(
      "Search for Term:\n Term: ",
      paste(term, collapse = ","),
      "\n\t Model: ",
      paste(formula(model), collapse = ", ")
    )
  }
  return(all(term %in% all_vars0(formula(model))))
}

SplitEffect <- function(finalformula,
                        fullModelFormula,
                        F.Model,
                        data,
                        depVariable,
                        mandatoryVar = "Genotype",
                        ci_levels = .95) {
  Allargs <- all_vars0(fullModelFormula)[!all_vars0(fullModelFormula) %in% c(depVariable, mandatoryVar)]
  if (is.null(Allargs)) {
    message0("Nothing to split on ...")
    return(NULL)
  }
  isCat <- !is.continuous(data[, Allargs, drop = FALSE])
  args <- Allargs[isCat]
  argsCon <- if (length(Allargs[!isCat]) > 0) {
    Allargs[!isCat]
  } else {
    # 	NULL
    1
  }
  largs <- length(args)
  dname <- names(data)
  l <- NULL
  names <- c()
  counter <- 1
  if (largs > 0) {
    for (i in seq_len(largs)) {
      argComb <- combn(args, i, simplify = TRUE)
      for (j in seq_len(ncol(argComb))) {
        arg <- as.vector(argComb[, j])
        if (arg %in% dname &&
          !is.numeric(data[, arg]) &&
          termInTheModel(model = fullModelFormula, term = arg)) {
          message0(
            counter,
            ". Split on ",
            paste0(arg, collapse = " & "),
            " ..."
          )
          newModel <- update(
            reformulate(
              response = depVariable,
              termlabels = c(
                paste(
                  mandatoryVar,
                  arg,
                  collapse = "*",
                  sep = "*"
                ),
                argsCon
              ),
              intercept = TRUE
            ),
            paste0(".~.-", mandatoryVar)
          )
          message0(
            "Checking the split model:\n\t",
            printformula(newModel)
          )

          l0 <- tryCatch(
            update(
              F.Model,
              newModel
            ),
            error = function(e) {
              message0(e, breakLine = FALSE)
              return(NULL)
            },
            warning = function(w) {
              message0(w, breakLine = FALSE)
              return(NULL)
            }
          )
          message0(
            "\tTested model: ",
            printformula(newModel),
            ifelse(!is.null(l0), " [Successful]", " [Failed]"),
            breakLine = FALSE
          )
          if (!is.null(l0)) {
            l0 <- intervalsCon(object = l0, lvls = ci_levels)
            l0$MainEffect <- arg
            l0$SplitFormula <- printformula(newModel)
            l[[counter]] <- l0
            names[counter] <- paste(
              paste0(mandatoryVar, collapse = "_"),
              paste0(arg, collapse = "."),
              collapse = "_",
              sep      = "_"
            )
            counter <- counter + 1
          }
        }
      }
    }
    if (length(names) > 0) {
      message0(
        "SplitEffects. Output names: ",
        paste0(names, collapse = ", ")
      )
      names(l) <- paste0(names)
    }
  }
  return(l)
}

dim0 <- function(...) {
  args <- list(...)
  r <- lapply(args, function(x) {
    if (is.null(dim(x))) {
      return(length(x))
    }
    dim(x)
  })
  unlist(r)
}

printformula <- function(formula, message = TRUE) {
  if (!is.null(formula)) {
    r <- paste01(format(formula, trim = TRUE, width = 0), collapse = "")
  } else {
    if (message) {
      message0("Ops! the formula is blank")
    }
    r <- NULL
  }
  return(r)
}
# Categorical effect size
# https://github.com/mpi2/stats_working_group/raw/master/PhenStatUserGuide/PhenStatUsersGuide.pdf p122
cat.eff.size <- function(xtb,
                         varName = NULL,
                         formula = NULL) {
  if (any(dim(xtb) < 1)) {
    r <- NULL
  } else {
    r <- max(apply(prop.table(xtb, margin = 2), 1, function(x) {
      max(dist(x, method = "maximum", diag = TRUE), na.rm = TRUE)
    }), na.rm = TRUE)
  }
  out <- list(
    value = r,
    variable = ifelse(!is.null(varName),
      varName,
      "Variable does not exist"
    ),
    model = printformula(RightFormula2LeftFormula(formula, removeResIfExists = TRUE)),
    type = "Proportion change",
    "percentage change" = NULL
  )
  return(out)
}

printVarFreqfromTable <- function(tbl) {
  if (is.null(tbl) || any(dim(tbl) < 1)) {
    return(NULL)
  }
  r <- pasteComma(paste0(names(tbl), "[", tbl, "]"))
  return(r)
}

GenderIncludedInAnalysis <- function(x, sexCol = "Sex") {
  r <- if (!sexCol %in% colnames(x)) {
    paste0(sexCol, " does not included in the input data")
  } else if (nlevels(x[, sexCol]) > 1) {
    paste0("Both sexes included; ", printVarFreqfromTable(table(x$Sex)))
  } else {
    paste0(
      "Only one sex included in the analysis; ",
      printVarFreqfromTable(table(x$Sex))
    )
  }
  return(r)
}

# Test engine
ctest <- function(x,
                  formula = NULL,
                  asset = NULL,
                  rep = 1500,
                  ci_levels = 0.95,
                  RRextraResults = NULL,
                  overallTableName = "Complete table",
                  InterLevelComparisions = TRUE,
                  ...) {
  xtb <- xtabs(
    formula = formula,
    data = x,
    drop.unused.levels = TRUE,
    na.action = "na.omit"
  )

  checkRes <- checkTableForFisherTest(
    xtb = xtb,
    asset = asset,
    check = c(1, 0)
  )
  checkResZeroVariation <- checkTableForFisherTest(
    xtb = xtb,
    asset = asset,
    check = c(0, 1)
  )
  if (!checkRes$passed) {
    r <- setNames(
      list(overall = list(
        p.value = NULL,
        effect = NULL
      )),
      overallTableName
    )
  } else {
    if (!checkResZeroVariation$passed || any(dim(xtb) < 2)) {
      r <- setNames(
        list(overall = list(
          p.value = 1,
          effect = NULL
        )),
        overallTableName
      )
    } else {
      r <- fisher.test1(
        x           = xtb,
        formula     = formula,
        ci_levels   = ci_levels,
        simulate.p.value = rep > 0,
        conf.int    = TRUE,
        conf.level  = ci_levels,
        B           = rep,
        overallTableName = overallTableName,
        InterLevelComparisions = InterLevelComparisions,
        ...
      )
    }
  }
  return(
    list(
      result = r,
      note = c(checkRes$message, checkResZeroVariation$message),
      table = xtb,
      input = x,
      formula = formula,
      RRextra = RRextraResults
    )
  )
}

# check table for fisher.test
checkTableForFisherTest <- function(xtb,
                                    asset = NULL,
                                    check = c(1, 0)) {
  message <- NULL
  if (!is.null(asset)) {
    if (asset <= dim(xtb)[3]) {
      xtb <- xtb[, , asset]
    } else {
      message <- c(message, "There are some empty levels in data")
      return(list(
        passed = FALSE,
        note = message
      ))
    }
  }
  if (check[1] &&
    (length(dim0(xtb)) < 2 ||
      (sum(margin.table(xtb, margin = 2) > 0) < 2 &&
        sum(margin.table(xtb, margin = 1) > 0) < 2))) {
    message <- c(
      message,
      "Contingency table with one level only or sum of margins less than 2"
    )
    return(list(
      passed = FALSE,
      note = message
    ))
  }

  if (check[2] &&
    sum(colSums(xtb) > 0) < 2) {
    message <- c(
      message,
      "Contingency table with one non-zero level"
    )
    return(list(
      passed = FALSE,
      note = message
    ))
  }

  return(list(
    passed = TRUE,
    note = message
  ))
}

# Bulletproof fisher.test
fisher.test0 <- function(x, formula, ci_levels, ...) {
  r <- tryCatch(
    do.call(fisher.test, listFun(
      list = list(x = x, ...),
      FUN = fisher.test
    )),
    error = function(e) {
      message0(e, breakLine = FALSE)
      return(NULL)
    },
    warning = function(w) {
      message0(w, breakLine = FALSE)
      return(NULL)
    }
  )
  r <- intervalsCat(r, ci_levels)
  r$effect <- cat.eff.size(x,
    varName = pasteComma(all_vars0(formula)[-1], truncate = FALSE),
    formula = formula
  )
  r$formula <- formula
  r$table <- x
  return(r)
}

# Fisher test with broken table
fisher.test1 <- function(x,
                         formula,
                         ci_levels,
                         overallTableName = "Complete table",
                         InterLevelComparisions = TRUE,
                         ...) {
  if (is.null(x)) {
    return(NULL)
  }

  nrx <- nrow(x)
  outList <- NULL
  if (InterLevelComparisions && nrx > 2) {
    message0(
      "\t\t testing sub tables in progress ...\n\t\t\t Total tests: ",
      ncombn(n = nrx, x = 2:nrx),
      "; ",
      pasteComma(names(attributes(x)$dimnames))
    )
    for (i in 2:nrx) {
      cbn <- combn(nrx, i)
      subtbl <- lapply(seq_len(ncol(cbn)), function(j) {
        r <- suppressMessages(fisher.test0(
          x = x[cbn[, j], ],
          formula = formula,
          ci_levels = ci_levels,
          ...
        ))
        if (nrow(cbn[, j, drop = FALSE]) == nrow(x)) {
          # this name is used in more places!
          r$data.name <- overallTableName
        } else {
          r$data.name <- paste(rownames(x[cbn[, j], ]), sep = ".", collapse = ".")
        }
        return(r)
      })

      outList <- c(outList, setNames(object = subtbl, lapply(subtbl, function(x) {
        x$data.name
      })))
    }
  } else {
    outList <- setNames(list(
      overall = fisher.test0(
        x         = x,
        formula   = formula,
        ci_levels = ci_levels,
        ...
      )
    ), overallTableName)
  }
  return(outList)
}

ncombn <- function(n, x, total = TRUE) {
  r <- lfactorial(n) - (lfactorial(x) + lfactorial(n - x))
  if (total) {
    return(sum(exp(r)))
  } else {
    return(exp(r))
  }
}


# Change with super extra care
# This is a complicated function for modeling all variation of the variables
AllTables <- function(dframe = NULL,
                      # dataframe
                      vars = NULL,
                      # list of categorical variables
                      cl = 0,
                      # lock the number of columns: works only if shrinke = TRUE
                      response.name = NULL,
                      # response name in the beginning of each variable name [only]
                      cols = NULL,
                      # filter on the columns of interest : works only if shrinke = TRUE
                      shrink = FALSE,
                      # Dichotomiing the final tables
                      Dichotomise = TRUE,
                      # This parameter is added for the MM! framework, adj = 0
                      adj = 1) {
  # remove no variation levels (columns)
  if (is.null(dframe)) {
    message0("Null data frame ")
    return(NULL)
  }

  cat <- vars
  lcat <- length(cat)
  # Make all tables
  l2 <- list()
  message0("\tSplitting in progress ...")
  for (i in unique(pmax(1, 1:(lcat - adj)))) {
    cb <- combn(cat, i)
    for (j in seq_len(ncol(cb))) {
      message0("\tSpliting on ", pasteComma(cb[, j], replaceNull = FALSE), " ...")
      out <- split(dframe,
        interaction(dframe[cb[, j]]),
        drop = TRUE
      )
      l2 <- c(l2, out)
    }
  }
  # Remove fixed value columns
  if (shrink) {
    message0("\tShrinking in progress ...")
    l2 <- lapply(l2, function(x) {
      # remove fixed value columns
      NonZeroFreq <- c(apply(x, 2, function(xx) {
        length(unique(na.omit(xx))) > 1
      }))
      NonZeroFreq [names(NonZeroFreq) %in% c("Freq", response.name)] <- TRUE
      r <- x[, NonZeroFreq, drop = FALSE]
      return(r)
    })
  }

  if (!is.null(cols)) {
    message0("\tKeeping variables of interest ...")
    l2 <- l2[as.logical(lapply(
      # keep certain columns
      l2,
      FUN = function(x) {
        all(cols %in% colnames(x))
      }
    ))]
  }

  # You do not want to get split on all combinations!
  if (all(cl > 0)) {
    # Fixed number of columns in each element of the list
    l2 <- l2[which(lapply(
      l2,
      FUN = function(x) {
        ncol(x) %in% cl
      }
    ) == TRUE)]
  }

  if (Dichotomise) {
    # Split the big tables into 2x2 tables
    message0("\tDichotomising the final tables ...")
    oblCols <- c(response.name, "Freq")
    l3 <- lapply(
      names(l2),
      FUN = function(z) {
        x <- l2[[z]]
        r <- list()
        nl <- names(x)[!names(x) %in% oblCols]
        if (length(nl) > 1) {
          cmbn <- combn(nl, length(nl) - 1)
          for (i in seq_len(ncol(cmbn))) {
            r[[i]] <- lapply(i, function(y) {
              x[, -which(names(x) %in% cmbn[, y]), drop = FALSE]
            })
            names(r[[i]]) <- paste(z, paste0(cmbn[, i], collapse = "...."), sep = "....")
          }
          return(unlist(r, recursive = FALSE))
        } else {
          return(x)
        }
      }
    )
    names(l3) <- names(l2)
    # Which sublists are sublevels?
    message0("\tFinalising the tables ....")
    k <- unlist(lapply(l3, function(x) {
      all(vapply(x, is.list, is.logical(1)))
    }), recursive = FALSE)
    f <- function(l) {
      names(l) <- NULL
      r <- unlist(l, recursive = FALSE, use.names = TRUE)
      return(r)
    }
    if (any(k)) {
      l3 <- c(l3[!k], f(l3[k]))
    }
  } else {
    l3 <- l2
  }
  return(l3)
}


FormulaHasIntercept <- function(formula) {
  if (is.null(formula)) {
    return(NULL)
  }
  attr(terms(as.formula(formula)), which = "intercept") > 0
}

FormulaHasResponse <- function(formula) {
  if (is.null(formula)) {
    return(NULL)
  }
  attr(terms(as.formula(formula)), which = "response") > 0
}
getResponseFromFormula <- function(formula) {
  if (is.null(formula)) {
    return(NULL)
  }
  if (attr(terms(as.formula(formula)), which = "response")) {
    all_vars0(formula)[1]
  } else {
    NULL
  }
}


formulaTerms <- function(formula,
                         response = FALSE,
                         intercept = FALSE) {
  if (!is.null(formula)) {
    r <- c(
      if (response && FormulaHasResponse(formula)) {
        getResponseFromFormula(formula)
      } else {
        NULL
      },
      attr(terms(as.formula(formula)), which = "term.labels"),
      if (intercept &&
        FormulaHasIntercept(formula)) {
        1
      } else {
        NULL
      }
    )
  }
  else {
    r <- NULL
  }
  return(r)
}

expand.formula <- function(formula) {
  reformulate(
    termlabels = labels(terms(formula)),
    response =
      if (attr(terms(formula), "response") > 0) {
        formula[[2]]
      } else {
        NULL
      }
  )
}

UnlistCall <- function(x) {
  as(unlist(x), "character")
}

ListOperation <- function(x, FUN = NULL) {
  cnames <- names(x)
  if (is.null(cnames)) {
    return(x)
  }
  x1 <- lapply(cnames, function(y) {
    ListOperation(x[[y]])
  })
  x1 <- FUN(x1)

  return(x1)
}

# You may want to use list.clean from rlist package
prunelist <- function(x) {
  message0("Pruning list in progress ...")
  r <- lapply(x, function(y) {
    if (is.null(y) || length(y) < 1) {
      NULL
    } else {
      y
    }
  })
  return(r)
}


ModelInReference <- function(model,
                             reference,
                             responseIncluded = FALSE,
                             veryLower = ~ Genotype + 1) {
  mo <- formulaTerms(formula = model, intercept = TRUE)
  re <- formulaTerms(formula = reference, intercept = TRUE)
  r <- re[re %in% mo]
  if (length(r) > 0) {
    out <- reformulate(
      termlabels = r,
      response = if (responseIncluded && FormulaHasResponse(model)) {
        all_vars0(model)[1]
      } else {
        NULL
      },
      intercept = TRUE
    )
    if (length(mo[!(mo %in% re)]) > 0) {
      message0(
        'Some terms in the "lower" model are ignored. See:\n\t',
        pasteComma(mo[!(mo %in% re)])
      )
      message0('The polished "lower": ', printformula(out))
    }
  } else {
    message0('An invalid "lower". It is set to: ', printformula(veryLower))
    out <- veryLower
  }
  return(out)
}

FERR_FullComparisionsMessage <- function(x) {
  message0("Optimisation level: ")
  message0("\tEstimation of all factor combination effects = ", x[1])
  message0("\tEstimation of inter level factors for the response = ", x[2])
}

TermInFormulaReturn <- function(formula,
                                term,
                                return,
                                active = TRUE,
                                not = NA,
                                debug = TRUE) {
  terms <- formulaTerms(formula = formula)
  if (debug) {
    message0(printformula(terms))
  }
  if (active && any(term %in% terms)) {
    return(return)
  } else {
    return(not)
  }
}

paste01 <- function(...) {
  r <- paste0(...)
  while (any(grepl(pattern = "  ", x = r))) {
    r <- gsub(
      pattern = "  ",
      x = r,
      replacement = " ",
      fixed = TRUE
    )
  }
  return(r)
}

greplM <- function(x = NULL, pattern = NULL, ...) {
  if (is.null(x)) {
    return(NULL)
  }
  if (is.null(pattern)) {
    return(x)
  }

  r <- rep(TRUE, length(x))
  for (p in pattern) {
    r <- r & grepl(pattern = p, x = x, ...)
  }
  return(r)
}

modelContrasts <- function(formula, data, ...) {
  if (is.null(formula) || is.null(data)) {
    return(NULL)
  }
  r <- colnames(model.matrix(formula, data, ...))
  return(r)
}

multiBatch <- function(data) {
  if ("Batch" %in% colnames(data)) {
    batchColumn <- na.omit(data[, "Batch"])
    if (length(levels(batchColumn)) > 1) {
      TRUE
    } else {
      FALSE
    }
  }
  else {
    FALSE
  }
}

unmatrix0 <- function(x, byrow = FALSE, sep = " ") {
  if (is.null(x)) {
    return(NULL)
  }

  rnames <- rownames(x)
  cnames <- colnames(x)
  if (is.null(rnames)) {
    rnames <- paste("r", seq_len(nrow(x)), sep = "")
  }
  if (is.null(cnames)) {
    cnames <- paste("c", seq_len(ncol(x)), sep = "")
  }
  nmat <- outer(rnames, cnames, paste, sep = sep)
  if (byrow) {
    vlist <- c(t(x))
    names(vlist) <- c(t(nmat))
  }
  else {
    vlist <- c(x)
    names(vlist) <- c(nmat)
  }
  return(vlist)
}

MoveResponseToRightOfTheFormula <- function(formula) {
  if (is.null(formula)) {
    message0("Ops! the formula is blank")
    return(NULL)
  }
  newFormula <- update.formula(
    old = formula,
    new = reformulate(
      response = NULL,
      termlabels = c(all_vars0(formula)[1], ".")
    )
  )
  out <- formula(delete.response(terms(newFormula)))
  message0("The input formula: ", printformula(formula))
  message0(
    "The reformatted formula for the algorithm: ",
    printformula(out)
  )
  return(out)
}

renameVariableNameInList <- function(list,
                                     name,
                                     replace = NULL,
                                     prefix,
                                     not = FALSE) {
  if (is.null(list)) {
    return(NULL)
  }

  if (is.null(replace)) {
    if (!not) {
      names(list)[names(list) %in% name] <- paste(prefix, names(list)[names(list) %in% name], sep = "_")
    } else {
      names(list)[!names(list) %in% name] <- paste(prefix, names(list)[!names(list) %in% name], sep = "_")
    }
  } else
  if (!not) {
    names(list)[names(list) %in% name] <- replace
  } else {
    names(list)[!names(list) %in% name] <- replace
  }

  return(list)
}

decimalplaces <- function(x) {
  if ((x %% 1) != 0) {
    nchar(strsplit(sub("0+$", "", as.character(x)), ".", fixed = TRUE)[[1]][[2]])
  } else {
    return(0)
  }
}

dataSignature <- function(formula, data, digits = 10) {
  a.vars <- all_vars0(formula)
  if (!is.null(formula) &&
    !is.null(data) &&
    !is.null(a.vars) &&
    nrow(data) > 1 &&
    sum(a.vars %in% names(data))) {
    v.vars <- a.vars[a.vars %in% names(data)]
    res0 <- lapply(v.vars, function(name) {
      d <- data[, name]
      if (is.numeric(d)) {
        paste0(
          name,
          ":[",
          pasteComma(
            paste0("n=", length(d)),
            paste0("mean=", round(mean(d, na.rm = TRUE), digits)),
            paste0("sd=", round(sd0(d, na.rm = TRUE), digits)),
            truncate      = FALSE,
            trailingSpace = FALSE
          ),
          "]"
        )
      } else {
        paste0(
          name,
          ":[",
          pasteComma(
            paste0(levels(d), "=", table(d)),
            truncate           = FALSE,
            trailingSpace      = FALSE
          ),
          "]"
        )
      }
    })

    overall <- mean(rowMeans(apply(data[, v.vars, drop = FALSE], 2, function(x) {
      r <- if (!is.numeric(x)) {
        as.integer(as.factor(x))
      } else {
        x
      }
      return(r)
    }), na.rm = TRUE))

    res0$overall <- paste0("Overall:[", overall, "]")
    res0$precision <- paste0("Precision:[", digits, "]")
    res <- pasteComma(
      sort(unlist(res0), decreasing = FALSE),
      truncate      = FALSE,
      trailingSpace = FALSE
    )
  } else {
    res <- paste0(
      "Something is wrong with the data/model. Make sure that the data is not null and the formula matches the data. [This is a random number ",
      randRegSeed(
        n = 1,
        decimal = TRUE,
        round = 10
      ),
      "]"
    )
  }
  return(res)
}

randRegSeed <- function(n = 1,
                        max = .Machine$double.xmax,
                        decimal = TRUE,
                        round = 10) {
  r <- runif(n, 1, as.numeric(Sys.time()) * 10000) %% max
  if (decimal) {
    r <- r - (r %% 1)
  }
  r <- round(r, digits = round)
  return(r)
}

as.numeric01 <- function(x) {
  if (!is.null(x)) {
    return(suppressWarnings(as.numeric(x)))
  } else {
    return(NULL)
  }
}

OpenStatsListLevels <- function(object, sep = "") {
  if (is.null(object)) {
    return(NULL)
  }
  l <- NULL
  SexLab <- "Sex"
  FemaleLab <- ifelse(
    is.na(object$input$OpenStatsList@dataset.values.female),
    "Female",
    object$input$OpenStatsList@dataset.values.female
  )
  MaleLab <- ifelse(
    is.na(object$input$OpenStatsList@dataset.values.male),
    "Male",
    object$input$OpenStatsList@dataset.values.male
  )
  SexLevels <- as.list(paste(SexLab, sort(c(
    FemaleLab, MaleLab
  ), decreasing = FALSE),
  sep = sep
  ))
  names(SexLevels) <- unlist(SexLevels)
  ##########
  GenotypeLab <- "Genotype"
  ControlLab <- ifelse(
    is.na(object$input$OpenStatsList@refGenotype),
    "control",
    object$input$OpenStatsList@refGenotype
  )
  MutantLab <- ifelse(
    is.na(object$input$OpenStatsList@testGenotype),
    "experimental",
    object$input$OpenStatsList@testGenotype
  )
  GenotypeLevels <- as.list(paste(GenotypeLab, sort(c(
    ControlLab, MutantLab
  ), decreasing = TRUE),
  sep = sep
  ))
  names(GenotypeLevels) <- unlist(GenotypeLevels)

  ##########
  BatchLab <- "Batch"
  ##########
  WeightLab <- "Weight"
  ##########
  l <- list(
    response = object$input$depVariable,
    LifeStage = list(
      LifeStage = "LifeStage",
      Early = "Early",
      Late = "Late",
      Levels = list(LifeStageEarly = "LifeStageEarly", LifeStageLate = "LifeStageLate")
    ),
    Sex = list(
      Sex = SexLab,
      Female = FemaleLab,
      Male = MaleLab,
      Levels = as.list(SexLevels)
    ),
    Genotype = list(
      Genotype = GenotypeLab,
      Control = ControlLab,
      Mutant = MutantLab,
      Levels = as.list(GenotypeLevels)
    ),
    Batch = BatchLab,
    Weight = WeightLab
  )
  return(l)
}

pasteCollon <- function(x) {
  r <- paste(x, collapse = ":", sep = ":")
  return(r)
}

# Modified from permutations in gtools package
permn <- function(n, r, v = seq_len(n), set = TRUE, repeats.allowed = FALSE) {
  if (mode(n) != "numeric" || length(n) != 1 || n < 1 ||
    (n %% 1) != 0) {
    stop("bad value of n")
  }
  if (mode(r) != "numeric" || length(r) != 1 || r < 1 ||
    (r %% 1) != 0) {
    stop("bad value of r")
  }
  if (!is.atomic(v) || length(v) < n) {
    stop("v is either non-atomic or too short")
  }
  if ((r > n) & repeats.allowed == FALSE) {
    stop("r > n and repeats.allowed=FALSE")
  }
  if (set) {
    v <- unique(sort(v))
    if (length(v) < n) {
      stop("too few different elements")
    }
  }
  # v0 <- vector(mode(v), 0)
  if (repeats.allowed) {
    sub <- function(n, r, v) {
      if (r == 1) {
        matrix(v, n, 1)
      } else if (n == 1) {
        matrix(v, 1, r)
      } else {
        inner <- Recall(n, r - 1, v)
        cbind(rep(v, rep(nrow(inner), n)), matrix(t(inner),
          ncol = ncol(inner), nrow = nrow(inner) * n,
          byrow = TRUE
        ))
      }
    }
  } else {
    sub <- function(n, r, v) {
      if (r == 1) {
        matrix(v, n, 1)
      } else if (n == 1) {
        matrix(v, 1, r)
      } else {
        X <- NULL
        for (i in seq_len(n)) {
          X <- rbind(X, cbind(v[i], Recall(n -
            1, r - 1, v[-i])))
        }
        X
      }
    }
  }
  sub(n, r, v[seq_len(n)])
}
CombineLevels <- function(...,
                          debug = TRUE,
                          len = 2) {
  x <- c(...)
  if (length(x) < 1) {
    return(x)
  }
  pm <- t(permn(n = length(x), r = len))
  r <- lapply(seq_len(ncol(pm)), function(i) {
    xx <- pm[, i]
    nxx <- length(xx)
    if (nxx < 1) {
      return(NULL)
    }
    r <- c(pasteCollon(x[xx]))
  })
  r <- as.vector(unlist(r))
  if (debug) {
    print(r)
  }
  return(r)
}

NullOrValue <- function(x, ReplaveValue = 10) {
  if (is.null(x)) {
    x <- ReplaveValue
  }
  return(x)
}

RRNewObjectAndFormula <- function(object,
                                  RRprop,
                                  formula,
                                  labels = NULL,
                                  depVarPrefix = NULL,
                                  refLevel = NULL,
                                  ### see right in the cut function
                                  right = TRUE) {
  allTerms <- all_vars0(formula)
  newobject <- object
  RRcutObject <- RRCut(
    object = object,
    prob = RRprop,
    depVariable = allTerms[1],
    labels = labels,
    depVarPrefix = depVarPrefix,
    right = right,
    refLevel = refLevel,
    lower = allTerms[2]
  )
  if (is.null(RRcutObject)) {
    return(NULL)
  }

  newobject <- RRcutObject$discObject
  newFormula <- replaceElementInFormula(
    formula = formula,
    pattern = allTerms[1],
    replace = RRcutObject$newdepVariable
  )
  return(
    list(
      newobject = newobject,
      newFormula = newFormula,
      newDepVariable = RRcutObject$newdepVariable,
      object = object,
      formula = formula,
      depVariable = allTerms[1],
      ###
      RRprop = RRprop,
      labels = labels,
      depVarPrefix = depVarPrefix,
      refLevel = RRcutObject$refLevel,
      empiricalQuantiles = RRcutObject$percentages
    )
  )
}

jitter0 <- function(x,
                    factor = 1,
                    amount = NULL,
                    upper = 1,
                    lower = .5,
                    maxtry = 1500) {
  for (i in seq_len(maxtry)) {
    if (i >= maxtry) {
      message0("\tNo solusion found for the specified RR_prop.")
    }
    xx <- jitter(
      x = x,
      amount = amount,
      factor = factor
    )
    if (min(xx, na.rm = TRUE) > lower && max(xx, na.rm = TRUE) < upper) {
      break
    }
  }
  return(xx)
}

addJitterToTheEntireData <- function(x, min = -1, max = 0) {
  if (is.null(x) || length(na.omit(x)) < 1) {
    return(x)
  }
  lx <- length(x)
  jitter <- sort(runif(n = lx, min = min, max = max), decreasing = FALSE)
  message0(
    "A small jitter (max value = ",
    min(jitter, na.rm = TRUE) -
      max(jitter, na.rm = TRUE),
    ") is added to the data"
  )
  o <- order(x)
  x <- x + jitter[o]
  return(x)
}

ExpandTails <- function(x,
                        amount,
                        upper = TRUE,
                        lower = TRUE) {
  xx <- na.omit(x)
  if (is.null(x) ||
    is.null(xx) ||
    length(xx) < 1) {
    return(x)
  }
  if (upper) {
    x[x == max(xx, na.rm = TRUE)] <- max(xx, na.rm = TRUE) + amount
  }
  if (lower) {
    x[x == min(xx, na.rm = TRUE)] <- min(xx, na.rm = TRUE) - amount
  }
  return(x)
}

ReplaceTails <- function(x,
                         lower = 0,
                         upper = 0,
                         add = FALSE) {
  xx <- na.omit(x)
  if (is.null(x) ||
    is.null(xx) ||
    length(xx) < 1) {
    return(x)
  }
  if (!add) {
    if (upper) {
      x[which.max(x)] <- upper[1]
    }
    if (lower) {
      x[which.min(x)] <- lower[1]
    }
  } else {
    x <- c(lower, x, upper)
  }
  return(sort(unname(x)))
}

CheckTheValidityOfTheRRLower <- function(lower, data, depvar, minLevels = 2) {
  if (is.null(data) || is.null(lower)) {
    stop("~> Null data or the `Reference variable`. Please check the input data or the `Reference variable`")
  }
  if (is.null(depvar) || !depvar %in% names(data)) {
    stop("~> dependent variable does not exist in the data")
  }
  if (!is.numeric(data[, depvar])) {
    stop("~> dependent variable must be numeric")
  }
  if (is.null(lower) || !lower %in% names(data)) {
    stop("~> `Reference variable` does not exist in the data")
  }
  if (is.numeric(data[, lower])) {
    stop("~> `Reference variable` must be a factor")
  }
  if (nlevels(as.factor(data[, lower])) < minLevels) {
    stop("~> `Reference variable` must have at least two levels")
  }
  return(TRUE)
}

RRGetTheLeveLWithMaxFrequency <- function(data, lower) {
  tbl <- table(data[, lower])
  maxTbl <- tbl %in% max(tbl, na.rm = TRUE)[1]
  r <- names(tbl[maxTbl])
  if (length(r) > 1) {
    message0(
      "\tMore than one variable with the highest frequency detected. See:\n\t Level(frequency): ",
      pasteComma(paste0(r, "(", tbl[maxTbl], ")"), truncate = FALSE),
      "\n\t\t The first one (`",
      r[1],
      "`) would be used."
    )
  } else {
    message0("\t\tNominated level(frequency) = ", pasteComma(paste0(r, "(", tbl[maxTbl], ")")))
  }
  return(r[1])
}

ExtractDatasetPLIfPossible <- function(x) {
  if (class(x) %in% c("PhenList", "OpenStatsList")) {
    x <- x@datasetPL
  }
  return(x)
}

RRCut <- function(object,
                  prob = .95,
                  depVariable = "data_point",
                  labels = NULL,
                  depVarPrefix = NULL,
                  right = TRUE,
                  refLevel = NULL,
                  lower = "Genotype",
                  decimal = 8) {
  data <- object
  if (class(data) %in% c("PhenList", "OpenStatsList")) {
    refLevel <- data@refGenotype
    data <- data@datasetPL
  }
  if (!CheckTheValidityOfTheRRLower(lower = lower, data = data, depvar = depVariable)) {
    return(NULL)
  }
  if (prob == .5) {
    message0("\t`prop` must be different from 0.5")
    return(NULL)
  }
  if (is.null(refLevel)) {
    message0("\tReference level left blank, then the dominate level will be set as the reference level")
    refLevel <- RRGetTheLeveLWithMaxFrequency(data = data, lower = lower)
  } else {
    message0("Reference level is set to `", refLevel, "`")
  }
  # Preparation ...
  JitterPrecision <- 4 + 1 * decimalplaces(min(data[, depVariable], na.rm = TRUE))
  data$data_point_discretised <- NA
  controls <- subset(
    data,
    data[, lower] %in% refLevel
  )
  mutants <- subset(
    data,
    !(data[, lower] %in% refLevel)
  )
  prb <- unique(c(0, prob, 1))
  qntl <- ReplaceTails(
    x = quantile(
      x = controls[, depVariable],
      probs = prb,
      na.rm = TRUE
    ),
    upper = max(data[, depVariable], na.rm = TRUE)[1] + 1,
    lower = min(data[, depVariable], na.rm = TRUE)[1] - 1,
    add = FALSE
  )
  if (sum(duplicated(qntl))) {
    message0(
      "\t* duplicates in quantiles detected, then small ",
      "(precision = 4 + minimum data precision) ",
      "jitters will be added to quantiles.\n\t\tQauntiles: ",
      pasteComma(round(qntl, 5), replaceNull = FALSE)
    )
    message0("\tJitter max precision (decimals) = ", JitterPrecision)
    qntl[duplicated(qntl)] <- jitter0(
      x = qntl[duplicated(qntl)],
      amount = 10^-JitterPrecision,
      upper = 1,
      lower = 0.5
    )
    qntl <- sort(qntl)
  }
  message0(
    "\tInitial quantiles for cutting the data \n",
    "\t\t\t Probs: ",
    pasteComma(round(prb, 3)),
    "\n\t\t\t N.reference: ",
    length(controls[, depVariable]),
    "\n\t\t\t Quantiles: ",
    pasteComma(round(qntl, 3), replaceNull = FALSE)
  )
  if (length(unique(qntl)) < 2) {
    message0("\tThe algorithm cannot specify non-unique quantiles.")
    return(NULL)
  }
  controls$data_point_discretised <- cut(
    x = controls[, depVariable],
    breaks = qntl,
    labels = if (is.null(labels)) {
      FALSE
    } else {
      labels
    },
    include.lowest = TRUE,
    right = right
  )
  ###
  tbc <- prop.table(table(controls$data_point_discretised))
  tbpc <- setNames(as.list(tbc), names(tbc))
  message0(
    "\t Detected percentiles in the data (",
    decimal,
    " decimals): ",
    pasteComma(paste(names(tbc), "=", round(tbc, 8)))
  )
  ###
  mutants$data_point_discretised <- cut(
    x = mutants [, depVariable],
    breaks = qntl,
    labels = if (is.null(labels)) {
      FALSE
    } else {
      labels
    },
    include.lowest = TRUE,
    right = right
  )
  newObj <- rbind(mutants, controls)
  newdepVariable <- pasteUnderscore(c(depVarPrefix, depVariable, "discretised"))
  names(newObj)[names(newObj) %in% "data_point_discretised"] <- newdepVariable
  newObj[, newdepVariable] <- as.factor(newObj[, newdepVariable])
  # message0('\tA new column is added to the data object: ', newdepVariable)
  return(
    list(
      object = object,
      discObject = newObj,
      newdepVariable = newdepVariable,
      depVariable = depVariable,
      percentages = tbpc,
      refLevel = refLevel,
      lower = lower
    )
  )
}

RRDiscretizedEngine <- function(data,
                                formula = data_point ~ Genotype + Sex + zygosity,
                                depVar = "data_point",
                                lower = "Genotype",
                                refLevel = NULL,
                                labels = c("Low", "NormalHigh"),
                                depVarPrefix = "Low",
                                right = TRUE,
                                prob = .95) {
  requireNamespace("rlist")
  l1 <- l2 <- l3 <- l4 <- NULL
  if (!CheckTheValidityOfTheRRLower(
    lower = lower,
    depvar = depVar,
    data = data
  )) {
    return(NULL)
  }
  vars <- names(data) %in% all.vars(formula)
  df <- data[, vars, drop = FALSE]
  df <- trimColsInDf(df = df)
  cat <- !is.continuous(df[, all.vars(formula), drop = FALSE])
  extra <- all.vars(formula)[cat &
    !all.vars(formula) %in% c(depVar, lower)]
  lextra <- length(extra)
  message0("Preparing the reference ranges ...")
  message0("Preparing the data for the variable: ", lower)

  l1 <- RRNewObjectAndFormula(
    object = data,
    RRprop = prob,
    formula = formula,
    labels = labels,
    depVarPrefix = depVarPrefix,
    refLevel = refLevel,
    right = right
  )

  if (lextra > 0) {
    for (i in seq_len(lextra)) {
      cbn <- combn(extra, i)
      for (j in seq_len(ncol(cbn))) {
        message0("\tSpliting on ", pasteComma(cbn[, j], replaceNull = FALSE), " ...")
        out <- split(df,
          interaction(df[cbn[, j]]),
          drop = TRUE
        )
        l2 <- c(l2, out)
      }
    }
  }
  if (!is.null(l2)) {
    l3 <- lapply(names(l2), function(name) {
      message0("Preparing the data for the combined effect: ", name)
      x <- l2[[name]]
      out <- RRNewObjectAndFormula(
        object = x,
        RRprop = prob,
        formula = formula,
        labels = labels,
        depVarPrefix = depVarPrefix,
        refLevel = refLevel,
        right = right
      )
    })
    names(l3) <- names(l2)
  }
  l4 <- c(list(l1), l3)
  #### must improve in future
  if (!is.null(l4)) {
    names(l4) <- ifelse(
      nzchar(names(l4)),
      paste(depVar, lower, names(l4), sep = "."),
      paste(depVar, lower, sep = ".")
    )
  }
  l4 <- list.clean(l4)
  return(l4)
}

RRextra <- function(object,
                    prob = .95,
                    depVariable = "data_point") {
  if (prob <= .5) {
    message0('"prob" must be greater than 0.5')
    return(NULL)
  }
  # Preparation ...
  controls <- subset(
    object@datasetPL,
    object@datasetPL$Genotype %in% object@refGenotype
  )
  mutants <- subset(
    object@datasetPL,
    object@datasetPL$Genotype %in% object@testGenotype
  )
  prb <- unique(c(0, 1 - prob, prob, 1))
  qntl <- quantile(x = controls[, depVariable], probs = prb)
  cutsC <- cut(x = controls[, depVariable], breaks = qntl)

  message0("Creating control cuts ...")
  XclassValue <- tapply(controls[, depVariable], cutsC, function(x) {
    (x)
  })
  message0("Creating mutant cuts ...")
  MclassValue <- lapply(
    XclassValue,
    FUN = function(xx) {
      sum(mutants$data_point %in% xx)
    }
  )
  CclassValue <- lapply(XclassValue, length)
  message0("Creating output tables ...")
  # Overall Table
  tbl <- rbind(unlist(CclassValue), unlist(MclassValue))
  dimnames(tbl) <- list(c("Control", "Mutant"), c("Low", "Normal", "High"))
  # Table Low
  tblLow <- cbind(tbl[, 1], rowSums(tbl[, 2:3]))
  dimnames(tblLow) <- list(c("Control", "Mutant"), c("Low", "Normal/High"))

  tblHigh <- cbind(tbl[, 1], rowSums(tbl[, seq_len(2)]))
  dimnames(tblHigh) <- list(c("Control", "Mutant"), c("Low/Normal", "High"))

  return(list(
    overall = tbl,
    tblLow = tblLow,
    tblHigh = tblHigh
  ))
}

TermInModelAndnLevels <- function(model,
                                  term = "LifeStage",
                                  data,
                                  threshold = 1) {
  if (is.null(model) ||
    is.null(data) || is.null(term) || !(term %in% colnames(data))) {
    return(FALSE)
  }
  # if (!is.null(model$correctd))
  # 	model      = model$correctd
  r <- termInTheModel(
    model = model,
    term = term,
    message = FALSE
  ) &&
    colLevelsSimple(data, all_vars0(model)[1]) > threshold
  return(r)
}
#######################
# From PhenStat
columnChecks0 <- function(dataset,
                          columnName,
                          dataPointsThreshold = 4) {
  presence <- TRUE
  numeric <- FALSE
  levelsCheck <- 0
  variabilityThreshold <- 10
  # Test: dependent variable presence
  if (!(columnName %in% colnames(dataset))) {
    presence <- FALSE
  }
  else {
    columnOfInterest <- na.omit(dataset[, c(columnName), drop = FALSE])

    if (all(is.continuous(columnOfInterest))) {
      numeric <- TRUE
    }

    dataPointsSummary <- columnLevels(dataset, columnName)

    NoCombinations <- dataPointsSummary[3]
    variabilityThreshold <- NoCombinations
    for (i in seq_len(NoCombinations)) {
      if (dataPointsSummary[3 + i] >= dataPointsThreshold) {
        levelsCheck <- levelsCheck + 1
      }
    }
  }

  values <-
    c(presence, numeric, (levelsCheck >= variabilityThreshold))

  return(values)
}

colLevelsSimple <- function(dataset, colName) {
  return(length(unique(dataset[, colName])))
}

columnLevels <- function(dataset, columnName) {
  columnOfInterest <- na.omit(dataset[, c(columnName)])


  values <- c(length(columnOfInterest))

  # Test for the data points quantity for Genotype/sex combinations
  Genotype_levels <- levels(factor(dataset$Genotype))
  Sex_levels <- levels(factor(dataset$Sex))
  values <- append(values, length(levels(factor(columnOfInterest))))

  values <-
    append(values, length(Genotype_levels) * length(Sex_levels))

  for (i in seq_along(Genotype_levels)) {
    GenotypeSubset <-
      subset(dataset, dataset$Genotype == Genotype_levels[i])
    for (j in seq_along(Sex_levels)) {
      GenotypeSexSubset <- subset(
        GenotypeSubset,
        GenotypeSubset$Sex == Sex_levels[j]
      )

      columnOfInterestSubset <-
        na.omit(GenotypeSexSubset[, c(columnName)])

      values <- append(values, length(columnOfInterestSubset))
    }
  }
  return(values)
}

# From  capitalize {Hmisc}
capitalise <- function(string) {
  capped <- grep("^[A-Z]", string, invert = TRUE)
  substr(string[capped], 1, 1) <- toupper(substr(
    string[capped],
    1, 1
  ))
  return(string)
}

sort0 <- function(x, ...) {
  if (length(x) > 0) {
    x <- sort(x = x, ...)
  }
  return(x)
}

message0 <- function(...,
                     breakLine = TRUE,
                     capitalise = TRUE,
                     appendLF = TRUE,
                     active = TRUE) {
  if (active) {
    x <- paste0(..., collapse = "")
    if (breakLine) {
      nmessage <- unlist(strsplit(x = x, split = "\n"))
    } else {
      nmessage <- x
    }
    if (capitalise) {
      nmessage <- capitalise(nmessage)
    }
    message(paste(Sys.time(), nmessage, sep = ". ", collapse = "\n"),
      appendLF = appendLF
    )
  }
}

warning0 <- function(...,
                     breakLine = TRUE,
                     capitalise = TRUE) {
  x <- paste0(..., collapse = "")
  if (breakLine) {
    nmessage <- unlist(strsplit(x = x, split = "\n"))
  } else {
    nmessage <- x
  }
  if (capitalise) {
    nmessage <- capitalise(nmessage)
  }
  warning(paste(Sys.time(), nmessage, sep = ". ", collapse = "\n"))
}


extractFisherSubTableResults <- function(x, what = "p.value") {
  r <- lapply0(x, function(y) {
    y[what]
  })
  return(as.list0(r))
}

lapply0 <- function(X, FUN, ...) {
  if (is.null(X)) {
    return(NULL)
  }
  r <- lapply(X = X, FUN = FUN, ...)
  return(r)
}

all_vars0 <- function(x, ...) {
  if (is.null(x)) {
    return(NULL)
  }
  fif <- all.vars(formula(x), ...)
  if (length(fif) > 0) {
    return(fif)
  } else {
    return(NULL)
  }
}

intervalsCon <- function(object, lvls, ...) {
  if (is.null(object)) {
    return(object)
  }

  message0(
    "\tComputing the confidence intervals at the level of ",
    pasteComma(lvls),
    " ..."
  )
  ci <- lapply(lvls, function(x) {
    citerms <- if (is(object, "lme")) {
      c("all", "fixed", "var-cov")
    } else if (is(object, "gls")) {
      c("all", "coef", "var-cov")
    } else {
      c("all")
    }
    for (citerm in citerms) {
      intv <- tryCatch(
        expr = intervals(
          object = object,
          level  = x,
          which  = citerm,
          ...
        ),
        error = function(e) {
          message0(
            "\t  ~> Error in estimating the confidence intervals for `",
            citerm,
            "` term(s)"
          )
          # message0('\t ~> ', e, breakLine = FALSE)
          return(NULL)
        },
        warning = function(w) {
          message0(
            "\t  ~> Error in estimating the confidence intervals for `",
            citerm,
            "` term(s)"
          )
          # message0('\t ~> ', w, breakLine = FALSE)
          return(NULL)
        }
      )
      if (!is.null(intv)) {
        message0("\t CI for `", citerm, "` term(s) successfully estimated")
        break
      } else if (!citerm %in% tail(citerms, 1)) {
        message0("\tAdjustment applied. Retrying ....")
      } else {
        message0("CI estimation failed.")
      }
    }
    if (is.null(intv)) {
      return(NULL)
    } else {
      return(list(intervals = intv, level = x))
    }
  })
  if (!is.null(ci)) {
    names(ci) <- paste("CI_", lvls, sep = "")
  }
  object$intervals <- ci
  return(object)
}



intervalsCat <- function(object, lvls = .95, ...) {
  if (is.null(object)) {
    return(object)
  }

  # message0('\t\tReformatting the confidence interval in the level of ',
  # 				 pasteComma(lvls),
  # 				 ' ...')
  c0 <- object$conf.int
  c2 <- NULL
  if (!is.null(c0)) {
    ci <- list(
      intervals = list(
        lower = c0[1],
        est. = as.vector(object$estimate),
        upper = c0[2]
      ),
      level = attr(c0, "conf.level")
    )
  } else {
    ci <- NULL
  }
  c2$intervals <- ci
  if (!is.null(ci)) {
    names(c2) <- paste("CI_", lvls, sep = "")
  }
  object$interval <- c2
  return(object)
}

expandDottedFormula <- function(formula, data) {
  if (is.null(formula) || is.null(data)) {
    message0("Null formula or data")
    return(formula)
  }

  newformula <- formula(terms.formula(x = formula, data = data))
  message0(
    "Extended formula (if (*.:) characters included):\n\t ",
    printformula(newformula)
  )
  return(newformula)
}

removeSingleLevelFactors <- function(formula, data) {
  cat <- all_vars0(formula)[!is.continuous(data[, all_vars0(formula), drop = FALSE])]
  if (length(cat)) {
    FactsThatMustBeRemoved <- cat[lapply(data[, cat, drop = FALSE], function(x) {
      length(unique(na.omit(x)))
    }) <= 1]
    if (length(FactsThatMustBeRemoved)) {
      message0(
        "The below terms from the model are removed because they only contain one level:\n\t ",
        pasteComma(
          FactsThatMustBeRemoved,
          replaceNull = FALSE,
          truncate = FALSE
        )
      )
      formula <- update.formula(
        formula,
        reformulate0(
          termlabels = c(".", FactsThatMustBeRemoved),
          response = NULL,
          intercept = TRUE,
          sep = "-"
        )
      )
    }
  }
  return(formula)
}

removeSingleValueContinuousVariables <- function(formula, data) {
  cont <- all_vars0(formula)[is.continuous(data[, all_vars0(formula), drop = FALSE])]
  if (length(cont)) {
    VarsThatMustBeRemoved <- cont[lapply(data[, cont, drop = FALSE], function(x) {
      length(unique(na.omit(x)))
    }) <= 1]
    if (length(VarsThatMustBeRemoved)) {
      message0(
        "The below terms from the model are removed because they do not have any variation:\n\t ",
        pasteComma(
          VarsThatMustBeRemoved,
          replaceNull = FALSE,
          truncate    = FALSE
        )
      )
      formula <- update.formula(
        formula,
        reformulate0(
          termlabels = c(".", VarsThatMustBeRemoved),
          response = NULL,
          intercept = TRUE,
          sep = "-"
        )
      )
    }
  }
  return(formula)
}


InteractionAndValue <- function(x, VName = "p-value") {
  r <- list()
  if (!is.null(x)) {
    r[VName] <- x
    r["Criteria"] <- TRUE
  } else {
    r$Criteria <- FALSE
  }
  return(r)
}


modelSummaryPvalueExtract <- function(x,
                                      variable = "Genotype",
                                      anova = TRUE,
                                      what = c("Pr(>|z|)", "Pr(>Chi)", "p-value"),
                                      debug = TRUE,
                                      ci_display = FALSE) {
  if (is.null(x)) {
    return(NULL)
  }
  if (anova) {
    if (any(class(x) %in% "glm")) {
      mOrg <- anova0(x, test = "LRT")
    } else {
      mOrg <- anova0(x, type = "marginal")
    }
  } else {
    mOrg <- coef(summary1(x))
  }
  mOrg <- as.data.frame(mOrg)
  if (debug) {
    print(mOrg)
  }
  mSum <- mOrg[, colnames(mOrg) %in% what, drop = FALSE]
  if (!is.null(variable) && any(variable %in% rownames(mOrg))) {
    mSumFiltered <- mSum[rownames(mOrg) %in% variable, , drop = FALSE]
  } else {
    mSumFiltered <- NULL # NA?
  }
  if (ci_display && !is.null(mSumFiltered)) {
    fixedEffInters <- CheckWhetherNamesExistInListExtractCI(
      x = x$intervals[[1]]$intervals,
      lnames = c("coef", "fixed"),
      variable = variable
    )
    return(list(
      "Value"      = as.vector(unlist(mSumFiltered)),
      "Confidence" =  Matrix2List(x = fixedEffInters),
      "Level"      = ifelse(is.null(fixedEffInters), NULL, x$intervals[[1]]$level)
    ))
  } else {
    return(as.vector(unlist(mSumFiltered)))
  }
}

CheckWhetherNamesExistInListExtractCI <- function(x, lnames, variable, minusCol = 2) {
  if (!is.null(x) && any(names(x) %in% lnames)) {
    Toplst <- x[names(x) %in% lnames][[1]]
    Toplst <- Toplst[rownames(Toplst) %in% variable, -minusCol,
      drop = FALSE
    ]
    return(Toplst)
  } else {
    return(NULL)
  }
}

CatEstimateAndCI <- function(object) {
  if (is.null(object)) {
    return(NULL)
  }
  v <- object$interval[[1]]
  if (!is.null(v)) {
    return(list(
      "Value" = v$intervals$est.,
      "Confidence" = list(
        "Lower" = v$intervals$lower,
        "Upper" = v$intervals$upper
      ),
      "Level" = v$level
    ))
  } else {
    return(list(
      "Value" = "Only available for 2x2 tables",
      "Confidence" = list(
        "Lower" = NULL,
        "Upper" = NULL
      ),
      "Level" = NULL
    ))
  }
}

NormaliseDataFrame <- function(data,
                               colnames = NULL) {
  message0("Normalising the data.frame in progress ...")

  if (!is.null(colnames)) {
    colnames <- colnames[colnames %in% names(data)]

    if (length(colnames) < 1) {
      message0('No variable name found in the data. Please check "colnames" ...')
      return(data)
    }
  } else {
    message0(
      "No variable selected for normalisation. All numerical variables will be normalised."
    )
    colnames <- names(data)
  }
  ######
  data  [, colnames] <- as.data.frame(lapply(
    data[, colnames, drop = FALSE],
    FUN = function(x) {
      if (is.numeric(x) && length(unique(x)) > 1) {
        sdx <- sd0(x, na.rm = TRUE)
        r <- (x - mean(x, na.rm = TRUE)) / ifelse(!is.na(sdx) &&
          sdx > 0, sdx, 1)
      } else {
        r <- x
      }
      return(r)
    }
  ))
  return(data)
}

getlmeObjectConfigs <- function(obj) {
  l <- list(
    "Formula" = printformula(if (class(obj) %in% "lme") {
      obj$call$fixed
    } else {
      obj$call$model
    }, message = FALSE),
    "Random effect" = printformula(obj$call$random, message = FALSE)
  )
  return(l)
}

extractLmeTerms <- function(object) {
  if (is.null(object) || !is.null(object$messages)) {
    return(NULL)
  }
  r <- list(
    "Initial model" = getlmeObjectConfigs(object$output$Initial.Model),
    "Final model"   = getlmeObjectConfigs(object$output$Final.Model)
  )
  return(r)
}

extractFERRTerms <- function(object) {
  if (is.null(object) || !is.null(object$messages)) {
    return(NULL)
  }
  lnitial <- final <- list(
    formula = printformula(object$input$formula, message = FALSE),
    random_effect = NULL
  )
  final$formula <- printformula(RightFormula2LeftFormula(object$extra$Cleanedformula),
    message = FALSE
  )
  out <- list(
    "Initial formula" = lnitial,
    "Final formula" = final
  )

  if (!is.null(object$input$RRprop)) {
    out$"RR quantile" <- unlist(object$input$RRprop)
  }

  return(out)
}

RightFormula2LeftFormula <- function(formula, removeResIfExists = FALSE) {
  r <- if (!is.null(formula) &&
    length(all_vars0(formula)) > 0) {
    if (removeResIfExists && FormulaHasResponse(formula)) {
      formula <- update(formula, NULL ~ .)
    }

    update(
      as.formula(formula),
      reformulate0(
        termlabels = c(".", all_vars0(formula)[1]),
        response = all_vars0(formula)[1],
        sep = "-"
      )
    )
  } else {
    NULL
  }
  return(r)
}

sd0 <- function(x, ...) {
  if (!is.numeric(x)) {
    return(NA)
  }
  r <- if (length(na.omit(x)) > 1) {
    sd(x, ...)
  } else {
    0
  }
  return(r)
}

SummaryStats <- function(x,
                         formula,
                         # label = 'raw_data_summary_statistics',
                         lower = FALSE,
                         drop = TRUE,
                         sep = "_",
                         removeSpecialChars = FALSE,
                         replace = "_") {
  r <- NULL
  formula <- checkModelTermsInData(
    formula = formula,
    data = x,
    responseIsTheFirst = TRUE
  )
  depVar <- all_vars0(formula)[1]
  if (is.null(depVar)) {
    message0("Null response! check the formula and the data")
    return(NULL)
  }
  # do not move me
  if (any(dim(x) == 0)) {
    return("empty dataset")
  }

  cat <- all_vars0(formula)[!is.continuous(x[, all_vars0(formula), drop = FALSE])]
  if (length(cat) > 0) {
    lvls <- interaction(x[, cat], sep = sep, drop = drop)
  } else {
    lvls <- rep(depVar, nrow(x))
  }

  isNumeric <- is.numeric(x[, depVar])
  summaryT <- as.list(tapply(x[, depVar], INDEX = lvls, function(xx) {
    if (isNumeric) {
      c <- ifelse(length(na.omit(xx)) > 0, length(na.omit(xx)), 0)
      m <- ifelse(length(na.omit(xx)) > 0, mean(xx, na.rm = TRUE), NA)
      sd <- ifelse(length(na.omit(xx)) > 0, sd0(xx, na.rm = TRUE), NA)
      normTest <- normality.test0(xx)
      r <- list(
        "Count" = c,
        "Mean" = m,
        "SD" = sd,
        "Normality test" = normTest
      )
    } else {
      c <- ifelse(length(na.omit(xx)) > 0, length(na.omit(xx)), 0)
      r <- list(
        "Count" = c,
        "Mean" = NULL,
        "SD" = NULL
      )
    }
    return(r)
  }, default = -999.991233210123))
  ##
  fTmp <- function(isNum) {
    if (isNum) {
      r <- list(
        "Count" = 0,
        "Mean" = NA,
        "SD" = NA
      )
    } else {
      r <- list("Count" = 0)
    }
    return(r)
  }
  summaryT[summaryT %in% c(-999.991233210123)] <- fTmp(isNum = isNumeric)

  if (lower) {
    nnames <- tolower(names(summaryT))
  } else {
    nnames <- names(summaryT)
  }

  if (removeSpecialChars) {
    nnames <- RemoveSpecialChars(nnames, replaceBy = replace)
  }

  # r = list(lbl = summaryT)
  # names(r) = label
  r <- summaryT
  return(r)
}

replaceElementInFormula <- function(formula,
                                    pattern,
                                    replace) {
  as.formula(gsub(
    pattern = pattern,
    replacement = replace,
    x = formula,
    fixed = TRUE
  ))
}

RemoveSpecialChars <- function(x,
                               what = "[^0-9A-Za-z]",
                               replaceBy = "_",
                               message = FALSE) {
  what <- gsub(
    pattern = "]",
    x = what,
    replacement = paste0(replaceBy, "]")
  )
  if (message) {
    message0(
      "pattern: ",
      what,
      "; replaced by ",
      replaceBy
    )
  }
  r <- gsub(what, replaceBy, x, ignore.case = TRUE)
  ###
  if (any(nchar(r) < 1)) {
    RN <- RandomRegardSeed(1, stringOutput = TRUE, round = 8)
    r[nchar(r) < 1] <- paste("no_name",
      RN,
      sep = "_",
      collapse = "_"
    )
  }
  return(r)
}

RandomRegardSeed <- function(n = 1,
                             max = 10^9,
                             decimal = TRUE,
                             round = 5,
                             stringOutput = FALSE,
                             what = "[^0-9A-Za-z]",
                             replaceBy = "_") {
  r <- runif(n, 1, as.numeric(Sys.time()) * 10000) %% max
  if (decimal) {
    r <- r %% 1
  }
  r <- round(r, digits = round)
  if (stringOutput) {
    r <- RemoveSpecialChars(as.character(r), what = what, replaceBy = replaceBy)
  }
  return(r)
}


####
rndProce <- function(procedure = NULL) {
  if (length(procedure) < 1 ||
    is.null(procedure) ||
    # !(procedure %in% lop()) ||
    length(procedure) > 1) {
    message0(
      "Error in inputing the procedure symbol. Current input: ",
      ifelse(is.null(procedure), "NULL", procedure),
      ". The default random effect, (1 |  Batch), is selected."
    )
    return(reformulate(" 1 |  Batch", response = NULL, intercept = TRUE))
  }
  ##############################################################
  if (procedure == "TYPICAL") {
    random <- reformulate(" 1 |  Batch", response = NULL, intercept = TRUE)
  } else if (procedure == "ABR") {
    random <- reformulate(" 1 |  Batch", response = NULL, intercept = TRUE)
  } else if (procedure == "ACS") {
    random <- reformulate(" 1 |  Batch", response = NULL, intercept = TRUE)
  } else if (procedure == "ALZ") {
    random <- reformulate(" 1 |  Batch", response = NULL, intercept = TRUE)
  } else if (procedure == "BLK") {
    random <- reformulate(" 1 |  Batch", response = NULL, intercept = TRUE)
  } else if (procedure == "BWT") {
    random <- reformulate(" 1 |  Batch", response = NULL, intercept = TRUE)
  } else if (procedure == "CAL") {
    random <- reformulate(" 1 |  Batch", response = NULL, intercept = TRUE)
  } else if (procedure == "CBC") {
    random <- reformulate(" 1 |  Batch", response = NULL, intercept = TRUE)
  } else if (procedure == "CHL") {
    random <- reformulate(" 1 |  Batch", response = NULL, intercept = TRUE)
  } else if (procedure == "CSD") {
    random <- reformulate(" 1 |  Batch", response = NULL, intercept = TRUE)
  } else if (procedure == "DXA") {
    random <- reformulate(" 1 |  Batch", response = NULL, intercept = TRUE)
  } else if (procedure == "ECG") {
    random <- reformulate(" 1 |  Batch", response = NULL, intercept = TRUE)
  } else if (procedure == "ECH") {
    random <- reformulate(" 1 |  Batch", response = NULL, intercept = TRUE)
  } else if (procedure == "ELZ") {
    random <- reformulate(" 1 |  Batch", response = NULL, intercept = TRUE)
  } else if (procedure == "EVL") {
    random <- reformulate(" 1 |  Batch", response = NULL, intercept = TRUE)
  } else if (procedure == "EVM") {
    random <- reformulate(" 1 |  Batch", response = NULL, intercept = TRUE)
  } else if (procedure == "EVO") {
    random <- reformulate(" 1 |  Batch", response = NULL, intercept = TRUE)
  } else if (procedure == "EVP") {
    random <- reformulate(" 1 |  Batch", response = NULL, intercept = TRUE)
  } else if (procedure == "EYE") {
    random <- reformulate(" 1 |  Batch", response = NULL, intercept = TRUE)
  } else if (procedure == "FER") {
    random <- reformulate(" 1 |  Batch", response = NULL, intercept = TRUE)
  } else if (procedure == "GEL") {
    random <- reformulate(" 1 |  Batch", response = NULL, intercept = TRUE)
  } else if (procedure == "GEM") {
    random <- reformulate(" 1 |  Batch", response = NULL, intercept = TRUE)
  } else if (procedure == "GEO") {
    random <- reformulate(" 1 |  Batch", response = NULL, intercept = TRUE)
  } else if (procedure == "GEP") {
    random <- reformulate(" 1 |  Batch", response = NULL, intercept = TRUE)
  } else if (procedure == "GPL") {
    random <- reformulate(" 1 |  Batch", response = NULL, intercept = TRUE)
  } else if (procedure == "GPM") {
    random <- reformulate(" 1 |  Batch", response = NULL, intercept = TRUE)
  } else if (procedure == "GPO") {
    random <- reformulate(" 1 |  Batch", response = NULL, intercept = TRUE)
  } else if (procedure == "GRS") {
    random <- reformulate(" 1 |  Batch", response = NULL, intercept = TRUE)
  } else if (procedure == "HEM") {
    random <- reformulate(" 1 |  Batch", response = NULL, intercept = TRUE)
  } else if (procedure == "HIS") {
    random <- reformulate(" 1 |  Batch", response = NULL, intercept = TRUE)
  } else if (procedure == "HWT") {
    random <- reformulate(" 1 |  Batch", response = NULL, intercept = TRUE)
  } else if (procedure == "IMM") {
    random <- reformulate(" 1 |  Batch", response = NULL, intercept = TRUE)
  } else if (procedure == "INS") {
    random <- reformulate(" 1 |  Batch", response = NULL, intercept = TRUE)
  } else if (procedure == "IPG") {
    random <- reformulate(" 1 |  Batch", response = NULL, intercept = TRUE)
  } else if (procedure == "OFD") {
    random <- reformulate(" 1 |  Batch", response = NULL, intercept = TRUE)
  } else if (procedure == "PAT") {
    random <- reformulate(" 1 |  Batch", response = NULL, intercept = TRUE)
  } else if (procedure == "VIA") {
    random <- reformulate(" 1 |  Batch", response = NULL, intercept = TRUE)
  } else if (procedure == "XRY") {
    random <- reformulate(" 1 |  Batch", response = NULL, intercept = TRUE)
  } else {
    random <- reformulate(" 1 |  Batch", response = NULL, intercept = TRUE)
  }
  return(random)
}

UniqueRatio <- function(x, ratio = TRUE) {
  if (length(x) < 1 ||
    length(na.omit(x)) < 1) {
    return(0)
  }
  x <- na.omit(x)
  n <- length(x)
  un <- length(unique(x))
  r <- ifelse(n > 0, un / n, 0)
  if (ratio) {
    r <- round(r * 100)
  }
  return(r)
}

mean0 <- function(x, na.rm = TRUE) {
  if (length(x) > 0 &&
    is.numeric(x)) {
    return(mean(x, na.rm = na.rm))
  } else {
    return(NULL)
  }
}

sd01 <- function(x, na.rm = TRUE) {
  if (length(x) > 0 &&
    is.numeric(x)) {
    return(sd0(x, na.rm = na.rm))
  } else {
    return(NULL)
  }
}

normality.test0 <- function(x, ...) {
  if (!is.null(x) &&
    is.numeric(x) &&
    length(x) > 3 &&
    !is.na(sd0(x, na.rm = TRUE)) &&
    length(unique(na.omit(x))) > 3 &&
    sd0(x, na.rm = TRUE) > 0) {
    #################### Shapiro
    if (length(x) < 5000) {
      r <- list(
        "P-value"          = shapiro.test(x, ...)$p.value,
        "Unique N"         = length(unique(na.omit(x))),
        "N"                = length(x),
        "Unique N/N percent" = UniqueRatio(x),
        "Mean"             = mean0(x),
        "SD"               = sd01(x),
        "Test"             = "Shapiro",
        "Note"             = "Cautions require when too many duplicates exist in data."
      )
    } else {
      #################### Kolmogorov-Smirnov
      precision <- 4 + decimalplaces(x = min(x, na.rm = TRUE))
      r <- list(
        "P-value" = ks.test(
          x = jitter(
            x = x,
            amount = precision
          ),
          y = "pnorm",
          alternative = "two.sided",
          ...
        )$p.value,
        "Unique N" = length(unique(na.omit(x))),
        "N" = length(x),
        "Unique N/N percent" = UniqueRatio(x),
        "Mean" = mean0(x),
        "SD" = sd01(x),
        "Test" = "Kolmogorov-Smirnov",
        "Note" = paste0(
          "Small jitter (precision = ",
          precision,
          " decimals) added to possible ties (duplicates)."
        )
      )
    }
  } else {
    r <- list(
      "P-value"            = NULL,
      "Unique N"           = length(unique(na.omit(x))),
      "N"                  = length(x),
      "Unique N/N percent" = UniqueRatio(x),
      "Mean"     = mean0(x),
      "SD"       = sd01(x),
      "Test"     = "No test applied to this data. Please check the data for possible QC issues",
      "Note"     = NULL
    )
  }
  return(r)
}


QuyalityTests <- function(object,
                          levels = c("Genotype", "Sex", "LifeStage"),
                          list = TRUE,
                          noDataLab = "No data",
                          sep = "_",
                          collapse = "_") {
  if (is.null(object) ||
    length(levels) < 1) {
    message0("\tSkipped . No model or the levels in the data for the quality test")
    return(NULL)
  }

  r <- resid(object)
  d <- getData(object)
  levels <- levels[levels %in% names(d)]
  levels <- levels[!is.continuous(d[, levels, drop = FALSE])]
  counter <- 1
  flst <- NULL
  if (length(levels) > 0) {
    for (i in seq_along(levels)) {
      cmb <- combn(x = levels, i)
      for (j in seq_len(ncol(cmb))) {
        result <- tapply(r, as.list(d[, cmb[, j], drop = FALSE]), function(x) {
          normality.test0(x)
        })
        if (!is.null(result)) {
          flst[[counter]] <- result
          names(flst)[counter] <- paste(cmb[, j], collapse = collapse, sep = sep)
          counter <- counter + 1
        }
      }
    }
    flst$Overall <- normality.test0(x = r)
    if (list) {
      flst <- as.list(lapply(flst, function(f) {
        r <- f
        if (length(dim(f)) > 1) {
          r <- as.matrix(as.data.frame(f))
          r <- unmatrix0(r)
        }
        r[is.na(r)] <- noDataLab
        r <- as.list(r)
        return(r)
      }))
    }
    return(flst)
  } else {
    return(NULL)
  }
}

as.list0 <- function(x, ...) {
  if (!is.null(x)) {
    return(as.list(x, ...))
  } else {
    return(NULL)
  }
}

AllEffSizes <- function(object, depVariable, effOfInd, data) {
  if (length(effOfInd) < 1 ||
    all(effOfInd == 1)) {
    message0("\tSkipped . No variable found for the effect size ...")
    return(NULL)
  }
  lst <- flst <- olst <- NULL
  data <- NormaliseDataFrame(data)
  effOfInd <- effOfInd[effOfInd %in% names(data)]
  counter1 <- counter2 <- 1
  if (!is.null(effOfInd) && length(effOfInd)) {
    cats <- effOfInd[!is.continuous(data[, effOfInd, drop = FALSE])]
    for (eff in effOfInd) {
      message0("\tLevel:", pasteComma(unlist(eff)))
      # 1. Main effect
      lstTmp <- eff.size(
        object = object,
        depVariable = depVariable,
        effOfInd = eff,
        errorReturn = NULL,
        data = data
      )
      if (!is.null(lstTmp)) {
        lst[[counter1]] <- lstTmp
        names(lst)[counter1] <- eff
        counter1 <- counter1 + 1
      }
      # 2. All subset interactions effect sizes
      if (length(cats) > 1) {
        for (j in seq_along(cats)) {
          cmbn <- combn(x = cats, m = j)
          for (k in seq_len(ncol(cmbn))) {
            interact <- interaction(data[, cmbn[, k], drop = FALSE], sep = " ", drop = TRUE)
            for (lvl in levels(interact)) {
              message0("\tLevel:", pasteComma(unlist(lvl)))
              olstTmp <- eff.size(
                object = object,
                data = droplevels0(data[interact %in% lvl, , drop = FALSE]),
                depVariable = depVariable,
                errorReturn = NULL,
                effOfInd = eff
              )
              if (!is.null(olstTmp)) {
                olst[[counter2]] <- olstTmp
                names(olst)[counter2] <- paste(eff, lvl, sep = "_")
                counter2 <- counter2 + 1
              }
            }
          }
        }
      }
    }
    flst <- as.list(c(lst, olst))
    return(flst)
  } else {
    return(NULL)
  }
}

extractAICc <- function(fit, scale = FALSE, k = NULL, ...) {
  requireNamespace("AICcmodavg")
  # k=0 is just loglikelihood or -2*logLik(fit)
  if (k == 0) {
    extractAIC(
      fit = fit,
      scale = scale,
      k = 0,
      ...
    )
  } else {
    edf <- AICc(mod = fit, return.K = TRUE, ...)
  }

  AIC <- AICc(mod = fit, return.K = FALSE, ...)
  l <- list(edf = edf, AIC = AIC)
  return(unlist(l))
}

stepAIC0 <- function(object,
                     scope,
                     scale = 0,
                     direction = c(
                       "both", "backward",
                       "forward"
                     ),
                     trace = 1,
                     keep = NULL,
                     steps = 1000,
                     use.start = FALSE,
                     k = 2,
                     ...) {
  if (is.null(object)) {
    message0("Null object in optimisation ...")
    return(object)
  }
  mydeviance <- function(x, ...) {
    dev <- deviance(x)
    if (!is.null(dev)) {
      dev
    } else {
      extractAICc(x, k = 0)[2L]
    }
  }
  cut.string <- function(string) {
    if (length(string) > 1L) {
      string[-1L] <- paste("\n", string[-1L], sep = "")
    }
    string
  }
  re.arrange <- function(keep) {
    namr <- names(k1 <- keep[[1L]])
    namc <- names(keep)
    nc <- length(keep)
    nr <- length(k1)
    array(unlist(keep, recursive = FALSE), c(nr, nc), list(
      namr,
      namc
    ))
  }
  step.results <- function(models, fit, object, usingCp = FALSE) {
    change <- sapply(models, "[[", "change")
    rd <- sapply(models, "[[", "deviance")
    dd <- c(NA, abs(diff(rd)))
    rdf <- sapply(models, "[[", "df.resid")
    ddf <- c(NA, abs(diff(rdf)))
    AIC <- sapply(models, "[[", "AIC")
    heading <- c(
      "Stepwise Model Path \nAnalysis of Deviance Table",
      "\nInitial Model:",
      deparse(formula(object)),
      "\nFinal Model:",
      deparse(formula(fit)),
      "\n"
    )
    aod <- if (usingCp) {
      data.frame(
        Step = change,
        Df = ddf,
        Deviance = dd,
        `Resid. Df` = rdf,
        `Resid. Dev` = rd,
        Cp = AIC,
        check.names = FALSE
      )
    } else {
      data.frame(
        Step = change,
        Df = ddf,
        Deviance = dd,
        `Resid. Df` = rdf,
        `Resid. Dev` = rd,
        AIC = AIC,
        check.names = FALSE
      )
    }
    attr(aod, "heading") <- heading
    class(aod) <- c("Anova", "data.frame")
    fit$anova <- aod
    fit
  }
  Terms <- terms(object)
  object$formula <- Terms
  if (inherits(object, "lme")) {
    object$call$fixed <- Terms
  } else if (inherits(object, "gls")) {
    object$call$model <- Terms
  } else {
    object$call$formula <- Terms
  }
  if (use.start) {
    warning("'use.start' cannot be used with R's version of 'glm'")
  }
  md <- missing(direction)
  direction <- match.arg(direction)
  backward <- direction == "both" | direction == "backward"
  forward <- direction == "both" | direction == "forward"
  if (missing(scope)) {
    fdrop <- numeric()
    fadd <- attr(Terms, "factors")
    if (md) {
      forward <- FALSE
    }
  }
  else {
    if (is.list(scope)) {
      fdrop <- if (!is.null(fdrop <- scope$lower)) {
        attr(terms(update.formula(object, fdrop)), "factors")
      } else {
        numeric()
      }
      fadd <- if (!is.null(fadd <- scope$upper)) {
        attr(terms(update.formula(object, fadd)), "factors")
      }
    }
    else {
      fadd <- if (!is.null(fadd <- scope)) {
        attr(terms(update.formula(object, scope)), "factors")
      }
      fdrop <- numeric()
    }
  }
  models <- vector("list", steps)
  if (!is.null(keep)) {
    keep.list <- vector("list", steps)
  }
  n <- nobs(object, use.fallback = TRUE)
  fit <- object
  bAIC <- extractAICc(fit, scale, k = k, ...)
  edf <- bAIC[1L]
  bAIC <- bAIC[2L]
  if (is.na(bAIC)) {
    stop("AIC is not defined for this model, so 'stepAIC' cannot proceed")
  }
  if (bAIC == -Inf) {
    stop("AIC is -infinity for this model, so 'stepAIC' cannot proceed")
  }
  nm <- 1
  Terms <- terms(fit)
  if (trace) {
    cat("Start:  AICc=",
      format(round(bAIC, 2)),
      "\n",
      cut.string(deparse(formula(fit))),
      "\n\n",
      sep = ""
    )
    utils::flush.console()
  }
  models[[nm]] <- list(
    deviance = mydeviance(fit),
    df.resid = n -
      edf,
    change = "",
    AIC = bAIC
  )
  if (!is.null(keep)) {
    keep.list[[nm]] <- keep(fit, bAIC)
  }
  usingCp <- FALSE
  while (steps > 0) {
    steps <- steps - 1
    AIC <- bAIC
    ffac <- attr(Terms, "factors")
    if (!is.null(sp <-
      attr(Terms, "specials")) && !is.null(st <- sp$strata)) {
      ffac <- ffac[-st, ]
    }
    scope <- factor.scope(ffac, list(add = fadd, drop = fdrop))
    aod <- NULL
    change <- NULL
    if (backward && length(scope$drop)) {
      aod <- dropterm0(
        fit,
        scope$drop,
        scale = scale,
        trace = max(
          0,
          trace - 1
        ),
        k = k,
        ...
      )
      rn <- row.names(aod)
      row.names(aod) <- c(rn[1L], paste("-", rn[-1L], sep = " "))
      if (any(aod$Df == 0, na.rm = TRUE)) {
        zdf <- aod$Df == 0 & !is.na(aod$Df)
        nc <- match(c("Cp", "AIC"), names(aod))
        nc <- nc[!is.na(nc)][1L]
        ch <- abs(aod[zdf, nc] - aod[1, nc]) > 0.01
        if (any(is.finite(ch) & ch)) {
          warning("0 df terms are changing AIC")
          zdf <- zdf[!ch]
        }
        if (length(zdf) > 0L) {
          change <- rev(rownames(aod)[zdf])[1L]
        }
      }
    }
    if (is.null(change)) {
      if (forward && length(scope$add)) {
        aodf <- addterm(
          fit,
          scope$add,
          scale = scale,
          trace = max(0, trace - 1),
          k = k,
          ...
        )
        rn <- row.names(aodf)
        row.names(aodf) <- c(rn[1L], paste("+", rn[-1L],
          sep = " "
        ))
        aod <- if (is.null(aod)) {
          aodf
        } else {
          rbind(aod, aodf[-1, , drop = FALSE])
        }
      }
      attr(aod, "heading") <- NULL
      if (is.null(aod) || ncol(aod) == 0) {
        break
      }
      nzdf <- if (!is.null(aod$Df)) {
        aod$Df != 0 | is.na(aod$Df)
      }
      aod <- aod[nzdf, ]
      if (is.null(aod) || ncol(aod) == 0) {
        break
      }
      nc <- match(c("Cp", "AIC"), names(aod))
      nc <- nc[!is.na(nc)][1L]
      o <- order(aod[, nc])
      if (trace) {
        print(aod[o, ])
        utils::flush.console()
      }
      if (o[1L] == 1) {
        break
      }
      change <- rownames(aod)[o[1L]]
    }
    usingCp <- match("Cp", names(aod), 0) > 0
    fit <- update(fit, paste("~ .", change), evaluate = FALSE)
    fit <- eval.parent(fit)
    nnew <- nobs(fit, use.fallback = TRUE)
    if (all(is.finite(c(n, nnew))) && nnew != n) {
      message0("\t Warning! number of rows in use has changed: remove missing values ...")
    }
    Terms <- terms(fit)
    bAIC <- extractAICc(fit, scale, k = k, ...)
    edf <- bAIC[1L]
    bAIC <- bAIC[2L]
    if (trace) {
      cat("\nStep:  AIC=",
        format(round(bAIC, 2)),
        "\n",
        cut.string(deparse(formula(fit))),
        "\n\n",
        sep = ""
      )
      utils::flush.console()
    }
    if (bAIC >= AIC + 1e-07) {
      break
    }
    nm <- nm + 1
    models[[nm]] <- list(
      deviance = mydeviance(fit),
      df.resid = n -
        edf,
      change = change,
      AIC = bAIC
    )
    if (!is.null(keep)) {
      keep.list[[nm]] <- keep(fit, bAIC)
    }
  }
  if (!is.null(keep)) {
    fit$keep <- re.arrange(keep.list[seq(nm)])
  }
  step.results(models = models[seq(nm)], fit, object, usingCp)
}

dropterm0 <-
  function(object,
           scope,
           scale = 0,
           test = c("none", "Chisq"),
           k = 2,
           sorted = FALSE,
           trace = FALSE,
           ...) {
    tl <- attr(terms(object), "term.labels")
    if (missing(scope)) {
      scope <- drop.scope(object)
    } else {
      if (!is.character(scope)) {
        scope <- attr(
          terms(update.formula(object, scope)),
          "term.labels"
        )
      }
      if (!all(match(scope, tl, 0L))) {
        stop("scope is not a subset of term labels")
      }
    }
    ns <- length(scope)
    ans <-
      matrix(
        nrow = ns + 1L,
        ncol = 2L,
        dimnames = list(c(
          "<none>",
          scope
        ), c("df", "AIC"))
      )
    ans[1, ] <- extractAICc(object, scale, k = k, ...)
    n0 <- nobs(object, use.fallback = TRUE)
    env <- environment(formula(object))
    for (i in seq_len(ns)) {
      tt <- scope[i]
      if (trace) {
        message(gettextf("trying - %s", tt), domain = NA)
        utils::flush.console()
      }
      nfit <- update(object, as.formula(paste("~ . -", tt)),
        evaluate = FALSE
      )
      nfit <- eval(nfit, envir = env)
      ans[i + 1, ] <- extractAICc(nfit, scale, k = k, ...)
      nnew <- nobs(nfit, use.fallback = TRUE)
      if (all(is.finite(c(n0, nnew))) && nnew != n0) {
        message0("\t Warning! number of rows in use has changed: remove missing values ...")
      }
    }
    dfs <- ans[1L, 1L] - ans[, 1L]
    dfs[1L] <- NA
    aod <- data.frame(Df = dfs, AIC = ans[, 2])
    o <- if (sorted) {
      order(aod$AIC)
    } else {
      seq_along(aod$AIC)
    }
    test <- match.arg(test)
    if (test == "Chisq") {
      dev <- ans[, 2L] - k * ans[, 1L]
      dev <- dev - dev[1L]
      dev[1L] <- NA
      nas <- !is.na(dev)
      P <- dev
      P[nas] <- safe_pchisq0(dev[nas], dfs[nas], lower.tail = FALSE)
      aod[, c("LRT", "Pr(Chi)")] <- list(dev, P)
    }
    aod <- aod[o, ]
    head <-
      c("Single term deletions", "\nModel:", deparse(formula(object)))
    if (scale > 0) {
      head <- c(head, paste("\nscale: ", format(scale), "\n"))
    }
    class(aod) <- c("anova", "data.frame")
    attr(aod, "heading") <- head
    aod
  }

safe_pchisq0 <- function(q, df, ...) {
  df[df <= 0] <- NA
  pchisq(q = q, df = df, ...)
}

lowHighList <- function(x, y, ...) {
  if (!is.null(x) || !is.null(y)) {
    r <- list("Low" = x, "High" = y, ...)
  } else {
    r <- list("Low" = x, "High" = y)
  }
  return(r)
}

factorise <- function(dataset, variables) {
  vars <- variables[variables %in% names(dataset)]
  if (length(vars) > 0) {
    for (v in vars) {
      dataset[, v] <- as.factor(dataset[, v])
    }
  }
  return(dataset)
}

OpenStatsListRelabling <- function(dataset, col, l1, rel1, l2, rel2) {
  if (col %in% names(dataset)) {
    dataset[, col] <- as.factor(dataset[, col])
    levels(dataset[, col]) <- capitalise(levels(dataset[, col]))
    if (!is.null(l1)) {
      levels(dataset[, col])[levels(dataset[, col]) %in% l1] <-
        rel1
    }
    if (!is.null(l2)) {
      levels(dataset[, col])[levels(dataset[, col]) %in% l2] <-
        rel2
    }
    lbls <- c(rel1, rel2)
    if (!is.null(lbls) && any(!levels(dataset[, col]) %in% lbls)) {
      message0(
        "There are some unused levels in `",
        col,
        "` that will be removed. \n\t Levels: ",
        pasteComma(levels(dataset[, col]))
      )
      RowIndFromCol <- dataset[, col] %in% lbls
      if (sum(RowIndFromCol) < 1) {
        message0(
          "\t  Preprocessing variable leads to an empty column \n\t   The variable renamed to `",
          paste0(col, "_labels`")
        )
        names(dataset)[names(dataset) %in% col] <- paste0(col, "_labels")
      } else {
        dataset <- droplevels(subset(dataset, RowIndFromCol))
      }
    }
  }
  return(dataset)
}

LevelsAsFacNumbered <- function(x, numbering = TRUE, sep = ". ") {
  if (is.null(x) ||
    length(x) < 1) {
    return(x)
  }
  r <- levels(as.factor(x))
  nr <- length(r)
  if (numbering) {
    r <- paste(seq_len(nr), r, sep = sep)
  }
  return(r)
}

numberingX <- function(x, sep = ". ") {
  if (length(x) < 1) {
    return(x)
  }
  lx <- length(x)
  r <- paste(seq_len(lx), x, sep = sep)
  return(r)
}

checkSummary <- function(dataset, var, numbering = TRUE, ...) {
  lvls <- NULL
  if (is.null(dataset) ||
    is.null(var) ||
    length(na.omit(dataset[, var])) < 1) {
    return(lvls)
  }

  MissingPercent <- round(sum(is.na(dataset[, var])) / length(dataset[, var]) * 100)
  MissingNote <- ifelse(
    length(MissingPercent) > 0 && MissingPercent > 50,
    paste0(MissingPercent, "% [Alert! Too many missings]"),
    paste0(MissingPercent, "%")
  )
  if (is.factor(dataset[, var]) || is.character(dataset[, var])) {
    lvls <- paste0(
      "\t Levels (Total levels = ",
      nlevels(as.factor(dataset[, var])),
      ", missings = ",
      MissingNote,
      "): \n\t  ",
      pasteComma(
        LevelsAsFacNumbered(dataset[, var], numbering = numbering),
        width = 250,
        sep = "\n\t  "
      )
    )
  } else {
    lvls <- paste0(
      "\t Summary:\n",
      "\t  mean      = ",
      mean(dataset[, var], na.rm = TRUE),
      "\n\t  sd        = ",
      sd0(dataset[, var], na.rm = TRUE),
      "\n\t  Missings  = ",
      MissingNote
    )
  }
  message0(lvls)
  return(invisible(lvls))
}

checkOpenStatsColumns <- function(dataset, vars) {
  allExist <- NULL
  if (length(vars)) {
    for (v in vars) {
      r <- v %in% names(dataset)
      message0(
        "Checking whether variable `",
        v,
        "` exists in the data ... \n\tResult = ",
        r
      )
      if (r) {
        checkSummary(dataset = dataset, var = v)
      }
      allExist <- c(allExist, r)
    }
  } else {
    allExist <- FALSE
  }
  return(allExist)
}


fIsBecome <- function(dataset,
                      is = "Assay.Date",
                      rename = "Batch",
                      ifexist = NULL) {
  if (is.null(is) || is.null(rename)) {
    return(dataset)
  }


  for (iss in is) {
    if (iss %in% colnames(dataset) &&
      ifelse(is.null(ifexist), TRUE, !iss %in% ifexist)) {
      ####
      if (!is.null(ifexist) &&
        ifexist %in% colnames(dataset) &&
        !ifexist %in% is) {
        colnames(dataset)[colnames(dataset) %in% rename] <- paste("Original", ifexist, sep = ".")
      }
      ###
      names(dataset)[names(dataset) %in% iss] <-
        rename
      message0(
        "Variable `",
        iss,
        "` renamed to `",
        rename, "`"
      )
      break
    }
  }
  return(dataset)
}

is.df.empty <- function(x) {
  if (is.null(x) || any(dim(x) < 1)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

CleanEmptyRecords <- function(x, vars) {
  if (is.df.empty(x)) {
    return(NULL)
  }
  vars1 <- vars[vars %in% names(x)]
  if (length(vars1)) {
    x <- x[all(!is.na(x[, vars1])) && all(x[, vars1] != ""), , drop = FALSE]
  }
  return(x)
}


fastBarnardextest <- function(x,
                              tail = 2,
                              prob = seq(10^-10, 1 - 10^-10, length.out = 101),
                              plot = FALSE) {
  if (is.null(x) ||
    !(is.table(x) ||
      is.matrix(x)) ||
    any(dim(x) != 2)) {
    stop("The input must be a 2 by 2 table/matrix")
  }

  fprob <- function(i, j, c1, c2) {
    n <- c1 + c2
    pa <- i / c1
    pb <- j / c2
    px <- (i + j) / n
    if (px == 0 || pa == pb) {
      return(0)
    } else {
      return((pa - pb) / sqrt(px * (1 - px) * ((1 / c1) + (1 / c2))))
    }
  }
  c <- colSums(x)
  r <- rowSums(x)
  c1 <- c[1]
  c2 <- c[2]
  n <- sum(c)
  pao <- x[1, 1] / c1
  pbo <- x[1, 2] / c2
  pxo <- r[1] / n
  TXO <- abs(pao - pbo) / sqrt(pxo * (1 - pxo) * (1 / c1 + 1 / c2))
  cbn <- matrix(c(rep(0:c1, each = c2 + 1), rep(0:c2, c1 + 1)), nrow = 2, byrow = TRUE)

  n1 <- lfactorial(c1)
  n2 <- lfactorial(c2)
  lprob <- log(prob)
  clprob <- log(1 - prob)
  Fact <- lfactorial(0:max(c1, c2, na.rm = TRUE)[1])
  ###################################################
  T <- sapply(seq_len(ncol(cbn)), function(col) {
    i <- cbn[1, col]
    j <- cbn[2, col]
    s <- n1 + n2 +
      (i + j) * lprob +
      (n - (i + j)) * clprob - sum(Fact[c(
        i + 1,
        j + 1,
        c1 - i + 1,
        c2 - j + 1
      )])
    t <- fprob(
      i = i,
      j = j,
      c1 = c1,
      c2 = c2
    )
    return(c(t = t, s = exp(s)))
  })
  r <- t(cbind(P = apply(T[-1, ], 1, function(x) {
    sum(x[T[1, ] >= TXO])
  }), prob = prob))
  ###################################################
  Nuisance.parameter <- unlist(r[2, ][which.max(r[1, ])])
  p.value <- unlist(r[1, ][which.max(r[1, ])])
  if (plot) {
    plot(
      r[2, ],
      r[1, ],
      type = "l",
      main = "Barnard's exact P-value",
      xlab = "Nuisance parameter",
      ylab = "P.value"
    )
    abline(v = Nuisance.parameter, col = 2, lwd = 2)
  }
  return(
    list(
      p.value = unname(min(tail * p.value, 1)),
      Nuisance.parameter = unname(Nuisance.parameter),
      Wald.Statistic = unname(TXO),
      tail = tail,
      seq = r
    )
  )
}

USerManualSlotName <- function(x, name = "OpenStats") {
  message("The structure of an ", name, ":")
  if (isS4(x)) {
    xn <- slotNames(x)
  } else {
    xn <- names(x)
  }

  if (length(xn) < 1) {
    return(NULL)
  }
  
  for (i in seq_along(xn)) {
    cat("    ", i, ". ", xn[i], "  \n")
  }
}


isVariableCategorical <- function(var = NULL, data = NULL) {
  if (is.null(data) ||
    is.null(var) ||
    !all(var %in% names(data))) {
    return(FALSE)
  }
  r <- ifelse(is.factor(data[, var]), TRUE, FALSE)
  return(r)
}

RemoveSexWithZeroDataPointInGenSexTableOnlyStatsPipelinenotExposed <- function(df = NULL,
                                                                               cols = c("Genotype", "Sex")) {
  if (is.null(df) ||
    nrow(df) < 1 ||
    !colExists(name = cols[1], data = df) ||
    !colExists(name = cols[2], data = df)) {
    return(df)
  }
  cols <- cols[cols %in% names(df)]
  if (length(cols) != 2) {
    message0("col parameter must have absolutely two values ...")
    return(df)
  }

  tbl <- table(df[, cols[1]], df[, cols[2]])
  if (length(dim(tbl)) > 1) {
    zc <- as.data.frame(tbl)
    zc <- zc[zc$Freq < 1, ]
    if (nrow(zc) > 0) {
      message0("Zero frequency values detected ...")
      names(zc)[seq_along(cols)] <- cols
      df <- df[!df[, cols[-1]] %in% zc[, cols[-1]], ]
      df <- droplevels(df)
    }
  }
  return(df)
}

MakeRRQuantileFromTheValue <- function(x, messages = TRUE) {
  if (length(x) < 1 ||
    length(na.omit(x)) < 1 ||
    !as.numeric(x) ||
    any(x > 1) ||
    any(x < 0)) {
    message0(
      "Quantile is not numeric [or not in [0,1] interval], value = `",
      pasteComma(x),
      "`"
    )
    return(x)
  }
  x.bck <- x
  x <- 1 - (1 - head(x, 1)) / 2
  if (messages) {
    message0(
      "The probability of the middle area in the distribution: ",
      pasteComma(x.bck)
    )
    message0(
      "\t Tails probability: ",
      min(x, 1 - x),
      "\n\t Formula to calculate the tail probabilities: 1-(1-x)/2, (1-x)/2 where x = ",
      x.bck
    )
  }
  return(max(x, 1 - x, na.rm = TRUE))
}

RRextraDetailsExtractor <- function(object,
                                    sep = " = ",
                                    collapse = ", ") {
  r <- NULL
  if (is.null(object) ||
    !is.null(object$messages)) {
    return(r)
  }

  r <- lapply(object$output$SplitModels, function(x) {
    lapply(x, function(y) {
      if (!is.null(y) && "RRextra" %in% names(y)) {
        y2 <- unlist(y$RRextra)
        if (is.null(y2)) {
          return(NULL)
        } else {
          return(paste(names(y2), y2, sep = sep, collapse = collapse))
        }
      }
    })
  })
  return(r)
}


trimColsInDf <- function(df, ...) {
  if (is.null(df) || !is.data.frame(df)) {
    return(df)
  }
  if (nrow(df) < 1) {
    return(df)
  }
  message0("Removing possible leading/trailing whitespace from the variables in the formula ...")
  df2 <- lapply(names(df), function(y) {
    x.bck <- x <- df[, y]
    if (is.numeric(x) || all(is.na(x)) || all(is.nan(x))) {
      x
    } else if (is.factor(x)) {
      levels(x) <- trimws(levels(x), ...)
    } else if (is.character(x)) {
      x <- trimws(x, ...)
    } else {
      x
    }
    if (!identical(x.bck, x)) {
      message0("\t Trim applied to variable `", y, "`")
    }
    return(x)
  })
  df2 <- droplevels(as.data.frame(df2))
  names(df2) <- names(df)
  return(df2)
}

FormatPvalueClassificationTag <- function(x, decimals = 4) {
  r <- ""
  if (is.null(x) || length(x) < 1) {
    return(r)
  }
  r <- if (min(x, na.rm = TRUE) >= 10^-decimals) {
    round(x, digits = decimals)
  } else {
    paste0(
      "< ",
      format(
        10^-decimals,
        digits = decimals,
        trim   = TRUE,
        nsmall = decimals,
        scientific = FALSE
      )
    )
  }
  return(r)
}

is.continuous <- function(x) {
  r <- vapply(x, is.numeric, is.logical(1))
  return(r)
}

totalListElements <- function(x) {
  if (is.null(x)) {
    return(0)
  }
  r <- sum(unlist(lapply(x, length)))
  return(r)
}

MainTitlePlusColon <- function(x = NULL) {
  x <- if (is.null(x)) {
    paste0(x, ": ")
  }
  return(x)
}

seq_along0 <- function(x, makeZero = TRUE) {
  s <- seq_along(x)
  if (makeZero && length(x) < 1) {
    s <- 1:0
  }
  return(s)
}
