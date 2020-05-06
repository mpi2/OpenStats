plot.NULL <- function(x, ...) {
  message0("No plot available for a NULL object")
}

plotFERR <- function(x, l1, l2, main, ...) {
  if (!is.null(x$output$SplitModels[[l1]][[l2]]$table)) {
    mosaicplot(x$output$SplitModels[[l1]][[l2]]$table, main = main, ...)
  } else {
    message0(
      "Can not find [",
      l1,
      " X ",
      l2,
      "] table"
    )
  }
}

###############################################
# Plot RR
###############################################
plot.OpenStatsRR <- function(x,
                             main = "Mosaic plot",
                             ask = FALSE,
                             mfrow = c(2, 2),
                             ...) {
  if (!is.null(x$messages) || is.null(x)) {
    message0("Due to error(s), no plot available")
    message0(x$messages)
    stop()
  }
  message0("Working on the plot ...")
  p <- par()
  par(ask = ask, mfrow = mfrow)

  Labels <- OpenStatsListLevels(x)
  LowRes <- pastedot("Low", Labels$response, "Genotype")
  HighRes <- pastedot("High", Labels$response, "Genotype")

  plotFERR(
    x = x,
    l1 = LowRes,
    l2 = Labels$Genotype$Genotype,
    main = main,
    ...
  )
  plotFERR(
    x = x,
    l1 = HighRes,
    l2 = Labels$Genotype$Genotype,
    main = main,
    ...
  )
  plotFERR(
    x = x,
    l1 = LowRes,
    l2 = Labels$Sex$Sex,
    main = main,
    ...
  )
  plotFERR(
    x = x,
    l1 = HighRes,
    l2 = Labels$LifeStage$LifeStage,
    main = main,
    ...
  )

  par(ask = p$ask, mfrow = p$mfrow)
}

###############################################
# Plot FE
###############################################
plot.OpenStatsFE <- function(x,
                             main = "Mosaic plot",
                             ask = FALSE,
                             mfrow = c(2, 2),
                             ...) {
  if (!is.null(x$messages) || is.null(x)) {
    message0("Due to error(s), no plot available")
    message0(x$messages)
    stop()
  }
  message0("Working on the plot ...")
  p <- par()
  par(ask = ask, mfrow = mfrow)

  Labels <- OpenStatsListLevels(x)
  plotFERR(
    x = x,
    l1 = Labels$response,
    l2 = Labels$Genotype$Genotype,
    main = main,
    ...
  )
  plotFERR(
    x = x,
    l1 = Labels$response,
    l2 = Labels$Sex$Sex,
    main = main,
    ...
  )
  plotFERR(
    x = x,
    l1 = Labels$Genotype$Genotype,
    l2 = Labels$Sex$Sex,
    main = main,
    ...
  )
  plotFERR(
    x = x,
    l1 = Labels$response,
    l2 = Labels$LifeStage$LifeStage,
    main = main,
    ...
  )
  par(ask = p$ask, mfrow = p$mfrow)
}

###############################################
# Plot MM
###############################################
plot.OpenStatsMM <- function(x,
                             main = "Final Model",
                             ask = FALSE,
                             mfrow = c(2, 2),
                             ...) {
  requireNamespace("car")
  if (!is.null(x$messages) || is.null(x)) {
    message0("Due to error(s), no plot available")
    message0(x$messages)
    stop()
  }
  message0("Working on the plot ...")
  fm <- x$output$Final.Model
  formula <- formula(fm)
  #transData <- applyFormulaToData(formula = formula, getData(fm))
  transData <- applyFormulaToData(formula = formula, x$input$OpenStatsList@datasetPL)
  n         <- nrow(transData$data)
  if (!identical(x = sortDataFrame(x$input$OpenStatsList@datasetPL),
                 y = sortDataFrame(getData(fm)))) {
    message0(
      "The input data is *NOT* the same as analysed dataset. ",
      "It is typically due to missing values in the input data"
    )
  }
  #
  p <- par()
  par(ask = ask, mfrow = mfrow)

  predR <- predict(fm)
  residR <- resid(fm)
  respShapiroTest  <- normality.test0(transData$data[, transData$names[1]],
                                      message=TRUE)
  residShapiroTest <- normality.test0(residR)
  plot(
    predR,
    residR,
    xlab = "Fitted values",
    ylab = "Residuals",
    main = main,
    ...
  )
  abline(h = 0, lwd = 3, lty = 2)
  densityPlot(
    na.omit(residR),
    xlab = ifelse(
      is.null(residShapiroTest$"P-value"),
      "Residuals",
      paste0(
        "Residuals - [",
        residShapiroTest$"Test",
        "] p-value = ",
        digit2Scientific(residShapiroTest$"P-value")
      )
    ),
    main = paste0(MainTitlePlusColon(main), "Density of the residuals"),
    rug = FALSE,
    method = "kernel",
    ...
  )
  qqPlot(
    na.omit(as.vector(residR)),
    ylab = "Residuals",
    main = paste0(MainTitlePlusColon(main), "Normal Q-Q of the residuals"),
    grid = FALSE,
    col.lines = 1,
    ...
  )
  # qqline(residR, ...)
  ptext <- paste0(transData$names[1], " (n = ", n, ")")
  densityPlot(
    transData$data[, transData$names[1]],
    main = "Density of the response",
    xlab = ifelse(
      is.null(respShapiroTest$"P-value"),
      ptext,
      paste0(
        ptext,
        " - [",
        respShapiroTest$"Test",
        "] p-value = ",
        digit2Scientific(respShapiroTest$"P-value")
      )
    ),
    rug = TRUE,
    ...
  )
  par(ask = p$ask, mfrow = p$mfrow)
}
