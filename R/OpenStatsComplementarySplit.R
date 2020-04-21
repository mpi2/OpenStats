OpenStatsComplementarySplit <- function(object = NULL,
                                        variables = c("Sex", "LifeStage"),
                                        debug = FALSE) {
  l <- NULL
  if (is.null(object)) {
    message0("~> Missing input object.")
    return(NULL)
  }
  if (!is.null(object$messages)) {
    message0("The input object is already failed. ")
    return(NULL)
  }
  vars <- unique(variables[variables %in% names(object$input$data)])
  nvars <- length(vars)
  if (!nvars) {
    message0(
      "Cannot find `variables` in the raw data.\n\t The input varaibles: ",
      pasteComma(variables)
    )
    return(NULL)
  }
  if (!object$input$method %in% "MM") {
    message0("Ineligible input object. The input object must be exported under MM framwork.")
    return(NULL)
  }

  message0("Split effects in progress ...")
  message0("Variable(s) to split:\n\t", pasteComma(vars))
  alTbls <- AllTables(
    dframe = as.data.frame(object$input$data),
    vars = vars,
    cl = 0,
    cols = NULL,
    response.name = NULL,
    shrink = FALSE,
    Dichotomise = 0,
    adj = 0
  )
  PLobj <- object$input$OpenStatsList
  ###
  l <- lapply(names(alTbls), function(x) {
    cat("\n")
    message0("Processing the levels: ", pasteComma(x))
    r <- tryCatch(
      expr = {
        PLobj@datasetPL <- droplevels(alTbls[[x]])
        r0 <- OpenStatsAnalysis(
          OpenStatsListObject = PLobj,
          method = object$input$method,
          MM_fixed = object$input$fixed,
          MM_random = object$input$random,
          MM_lower = object$input$lower,
          MM_weight = object$input$weight,
          MM_direction = object$input$direction,
          MM_checks = object$input$checks,
          MM_optimise = object$input$optimise,
          MMFERR_conf.level = object$input$ci_level,
          MM_BodyWeightIncluded = NULL,
          debug = debug
        )
      },
      warning = function(war) {
        message0("Submodel failed ...")
        message0(war)
        return(NULL)
      },
      error = function(err) {
        message0("Submodel failed ...")
        message0(err)
        return(NULL)
      }
    )
    ###
    if (!is.null(r) ||
      is.null(r$messages)) {
      message0("[Successful]")
    } else {
      message0("[Failed]")
    }

    return(r)
  })

  if (!is.null(l)) {
    names(l) <- names(alTbls)
  }

  class(l) <- "OpenStatsComplementarySplit"
  return(invisible(l))
}


plot.OpenStatsComplementarySplit <- function(x,
                                             main = "Final Model",
                                             ask = FALSE,
                                             mfrow = c(2, 2),
                                             ...) {
  R <- NULL
  if (is.null(x)) {
    message0("No plot available for a NULL object")
    return(R)
  }
  ObjectNames <- names(x)
  if (length(ObjectNames) > 0) {
    lapply(ObjectNames, function(xp) {
      plot(
        x[[xp]],
        main = paste0(xp, "\n", main),
        ask = ask,
        mfrow = mfrow,
        ...
      )
    })
  }
  return(invisible(R))
}

summary.OpenStatsComplementarySplit <- function(object, format = "rst", ...) {
  R <- NULL
  if (is.null(object)) {
    message0("No summary available for a NULL object")
    return(R)
  }
  ObjectNames <- names(object)
  if (length(ObjectNames) > 0) {
    R <- lapply(ObjectNames, function(xp) {
      message0("Summary for ", xp)
      summary(object[[xp]],
        format = format,
        caption = xp,
        ...
      )
    })
  }
  return(invisible(R))
}

print.OpenStatsComplementarySplit <- function(x, format = "rst", ...) {
  r <- summary.OpenStatsComplementarySplit(object = x, format = format, ...)
  return(invisible(r))
}
