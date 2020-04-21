# Categorical OPT core
crunner <- function(object,
                    formula = ~ category + Genotype + Sex + LifeStage,
                    rep = 1500,
                    method = NULL,
                    fullComparisions = TRUE,
                    noteToFinish = NULL,
                    ci_levels = .95,
                    RRextraResults = NULL,
                    overallTableName = "Complete table",
                    InterLevelComparisions = TRUE,
                    trimWC = TRUE,
                    ...) {
  requireNamespace("rlist")
  if (sum(method != c("FE", "RR")) == 0 ||
    is.null(object) ||
    is.null(all_vars0(formula)) ||
    length(all_vars0(formula)) < 2) {
    message0(
      "Improper method (",
      method,
      ") for the type of data, or the `formula/data` is not specified properly/left blank.\n\t(right sided) Formula: ",
      printformula(formula)
    )
    return(NULL)
  }
  # message0('FE framework in progress ...')
  sta.time <- Sys.time()
  message0("\tTop framework: ", method)
  message0(
    "Fisher exact test with ",
    ifelse(rep > 0, rep, "No"),
    " iteration(s) in progress ..."
  )

  data <- ExtractDatasetPLIfPossible(object)
  newFormula <- checkModelTermsInData(
    formula = formula,
    data = data,
    responseIsTheFirst = TRUE
  )
  missingInVariable(
    fixed = newFormula,
    data = data,
    threshold = 50
  )
  ####
  allTerms <- all_vars0(newFormula)
  if (length(allTerms) < 2) {
    message0(
      "No `response` to test or the right hand side of the `formula` not specified properly/left blank.\n\t(right sided) Formula: ",
      printformula(newFormula)
    )
    return(NULL)
  }
  ####
  if (trimWC) {
    data <- trimColsInDf(data[, allTerms, drop = FALSE])
  }
  lComplete <- NULL
  for (indx in seq_len(ifelse(fullComparisions, pmax(1, length(allTerms) - 1), 1))) {
    depVariable <- allTerms[indx]
    vars <- allTerms[-c(seq_len(indx))]
    l <- lcomb <- names <- alTbls <- NULL
    CmbiVars <- (length(vars) > 1)
    ####
    message0("Step ", indx, '. Testing "', depVariable, '"')
    ####
    newObject <- data
    Obj <- CheckMissing(
      newObject,
      reformulate(termlabels = depVariable, response = NULL)
    )
    newObject <- Obj$new.data
    ####
    if (!isVariableCategorical(var = depVariable, data = data)) {
      message0("\tThe response variable `", depVariable, "` is not categorical.")
      next
    }
    ####
    counter <- 1
    for (j in seq_along(vars)) {
      message0(
        "\tTesting for the main effect: ",
        pasteComma(vars[j], replaceNull = FALSE)
      )
      lt <- ctest(
        x = newObject,
        formula = reformulate(
          termlabels = c(depVariable, vars[j]),
          response = NULL
        ),
        rep = rep,
        ci_levels = ci_levels,
        RRextraResults = RRextraResults,
        overallTableName = overallTableName,
        InterLevelComparisions = InterLevelComparisions,
        ...
      )

      if (!is.null(lt) && !is.null(lt$result)) {
        l[[counter]] <- lt
        names <- c(names, vars[j])
        names(l) <- names
        counter <- counter + 1
      }
    }
    if (fullComparisions && CmbiVars) {
      message0("Combined effects in progress ...")
      alTbls <- AllTables(
        dframe = as.data.frame(xtabs(
          formula = newFormula,
          data = newObject
        )),
        vars = vars,
        cl = 3:(length(vars) + 2),
        cols = NULL,
        response.name = depVariable,
        shrink = TRUE
      )
      message0("Testing for the combined effects ... ")
      lcomb <- c(lcomb, lapply(alTbls, function(x) {
        ctest(
          x = x,
          formula = reformulate0(
            response = "Freq",
            termlabels = colnames(x)[!colnames(x) %in% "Freq"]
          ),
          rep = rep,
          ci_levels = ci_levels,
          RRextraResults = RRextraResults,
          overallTableName = overallTableName,
          InterLevelComparisions = InterLevelComparisions,
          ...
        )
      }))
      if (!is.null(lcomb)) {
        names(lcomb) <- paste(lapply(
          names(lcomb),
          FUN = function(nam) {
            # n1 = names(dimnames(lcomb[[nam]]$result[[overallTableName]]$table))
            n1 <- all_vars0(lcomb[[nam]]$formula)[-1]
            n2 <- paste0(n1[!n1 %in% depVariable], collapse = ".")
            return(n2)
          }
        ),
        names(lcomb),
        sep = "_"
        )
        l <- c(l, lcomb)
      }
    }
    # fix naming issue (make them short)
    if (!is.null(l) && length(l) > 0) {
      names(l) <- gsub("\\.{4}.*", "", names(l))
    }
    lComplete[[depVariable]] <- l
  }
  message0(
    "Total tested categories = ",
    length(lComplete),
    ": ",
    pasteComma(names(lComplete), replaceNullby = "-")
  )
  message0("\tTotal tests =  ", totalListElements(lComplete))
  message0(
    "FE framework ",
    noteToFinish,
    " executed in ",
    round(difftime(Sys.time(), sta.time, units = "sec"), 2),
    " second(s)."
  )
  OutR <- list(
    output = list(SplitModels = list.clean(lComplete)),
    input = list(
      OpenStatsList = object,
      data = data,
      depVariable = allTerms[1],
      rep = rep,
      method = method,
      formula = formula,
      ci_level = ci_levels,
      full_comparisions = c(fullComparisions, InterLevelComparisions)
    ),
    extra = list(
      missings = Obj$missings,
      UsedData = Obj$new.data,
      AllTable = alTbls,
      Cleanedformula = newFormula
    )
  )
  class(OutR) <- "OpenStatsFE"
  return(OutR)
}
