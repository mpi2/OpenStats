# Continuous OPT core
M.opt <- function(object = NULL,
                  fixed = NULL,
                  random = NULL,
                  lower = ~ Genotype + 1,
                  direction = "both",
                  trace = TRUE,
                  weight,
                  checks = c(1, 1, 1, 1),
                  method = "MM",
                  optimise = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
                  ci_levels = .95,
                  ...) {
  requireNamespace("nlme")
  if (!method %in% c("MM") ||
    is.null(all_vars0(fixed)) ||
    is.null(object)) {
    #####
    message0(
      "Improper method (",
      method,
      ") for the type of data, or the `formula/data` is left blank"
    )
    return(NULL)
  }
  # Future devs
  family <- poisson(link = "log")
  categorical <- FALSE
  # Start from here
  sta.time <- Sys.time()
  lCont <- lmeControl(
    opt = "optim",
    maxIter = 1500,
    msMaxIter = 1500
  )
  gCont <- glsControl(
    opt = "optim",
    singular.ok = TRUE,
    maxIter = 1500,
    msMaxIter = 1500
  )
  glCont <- glm.control(epsilon = 10^-36, maxit = 1500)
  G.Model <- FV.Model <- I.Model <- SplitModels <- EffectSizes <- F.Model <- OutR <- ResidualNormalityTest <- NULL
  VarHomo <- TRUE
  data <- RemoveDuplicatedColumnsFromDfandTrimWhiteSpace(x = object@datasetPL, formula = fixed)
  n <- nrow(data)
  fixed <- ModelChecks(
    fixed = fixed,
    data = data,
    checks = checks
  )
  CheckedRandom <- RandomEffectCheck(
    formula = random,
    data = data
  )
  allVars <- all_vars0(fixed)
  LifeStage <- "LifeStage" %in% allVars
  ##################
  Batch_exist <- !categorical && !is.null(CheckedRandom)
  initialFixed <- fixed
  Imdl <- ifelse(Batch_exist, "lme", ifelse(categorical, "glm", "gls"))
  ##################
  for (mdl in unique(c(Imdl, "gls"))) {
    message0(mdl, ": Fitting the full model ... ")
    fixed <- initialFixed
    fixedTerms <- formulaTerms(initialFixed)
    for (i in seq_along0(fixedTerms)) {
      I.Model <- tryCatch(
        expr = do.call(
          mdl,
          listFun(
            list = list(
              model = if (mdl == "glm") {
                TRUE
              } else {
                fixed
              },
              fixed = fixed,
              formula = fixed,
              family = family,
              random = CheckedRandom,
              data = data,
              na.action = na.omit,
              method = ifelse(mdl == "glm", "glm.fit", "REML"),
              weights = if (mdl != "glm") {
                weight
              } else {
                NULL
              },
              control = if (Batch_exist) {
                lCont
              } else {
                if (categorical) {
                  glCont
                } else {
                  gCont
                }
              },
              ...
            ),
            FUN = ifelse(Batch_exist, "lme", ifelse(categorical, "glm", "gls")),
            debug = TRUE
          )
        ),
        warning = function(war) {
          message0("* The full model failed with the warning (see below): ")
          message0("\t", war, breakLine = FALSE)
          return(NULL)
        },
        error = function(err) {
          message0("* The full model failed with the error (see below): ")
          message0("\t", err, breakLine = FALSE)
          return(NULL)
        }
      )
      ###########
      if (is.null(I.Model)) {
        removedTerm <- fixedTerms[length(fixedTerms) - (i - 1)]
        ltemp <- ModelInReference(model = lower, reference = fixed)
        if (!is.null(lower) &&
          !is.null(fixed) &&
          !is.null(ltemp) &&
          removedTerm %in% formulaTerms(ltemp)) {
          message0(
            "The following term did not removed: ",
            removedTerm,
            "\n\t the terms below will not remove from the model:\n\t ",
            printformula(formulaTerms(ltemp))
          )
          next
        }
        fixed <- update.formula(
          old = initialFixed,
          new = reformulate0(
            response = ".",
            termlabels = c(".", removedTerm),
            sep = "-"
          )
        )
        message0(
          "Round ",
          i,
          " of fixing the error. The followeing term will be removed: ",
          removedTerm,
          "\n\tNew formula: ",
          printformula(fixed)
        )
      } else {
        message0("\tThe full model successfully applied.")
        break
      }
    }
    if (!is.null(I.Model)) {
      break
    }
    if (Batch_exist && is.null(I.Model) && mdl %in% "lme") {
      message0(mdl, " failed. Retrying with a different model ...")
      Batch_exist <- FALSE
    }
  }
  ###########
  if (is.null(I.Model)) {
    message0("Full model failed ...")
  } else {
    I.Model <- intervalsCon(object = I.Model, lvls = ci_levels)
    message0(
      'The specified "lower" model: \n\t',
      ifelse(!is.null(lower), printformula(lower), "Null lower")
    )
    lowerCorrected <- ModelInReference(model = lower, reference = fixed)
    optimiseMessage(head(optimise, 3))
  }
  ###########
  if (!is.null(I.Model) && optimise[1] && !is.null(lowerCorrected)) {
    message0("\tThe direction of the optimisation (backward, forward, both): ", direction)
    message0("\tOptimising the model ... ")
    F.Model <- tryCatch(
      # stepAIC must have ML fit as the input
      # BIC/AIC is replaced with AICc then this parameter (k) does not work
      expr = stepAIC0(
        REML2ML(I.Model),
        trace     = trace,
        direction = direction,
        scope     = list(lower = lowerCorrected),
        na.action = na.omit,
        k = log(n)
      ),
      warning = function(war) {
        message0("\t * The optimisation failed with the warning (see below): ")
        message0("\t", war, breakLine = FALSE)
        return(NULL)
      },
      error = function(err) {
        message0("\t * The optimisation failed with the error (see below): ")
        message0("\t", err, breakLine = FALSE)
        return(NULL)
      }
    )
    ###########
    if (is.null(F.Model)) {
      message0("\tOptimisation did not apply")
      F.Model <- I.Model
      optimise[1] <- FALSE
    } else {
      message0("\tOptimised model: ", printformula(formula(F.Model)))
      F.Model <- ML2REML(F.Model)
      F.Model <- intervalsCon(object = F.Model, lvls = ci_levels)
      optimise[1] <- TRUE
    }
  } else {
    F.Model <- I.Model
    optimise[1] <- FALSE
  }
  ###########
  if (!is.null(F.Model)) {
    if (optimise[2] && !is.null(weight) && !(mdl %in% "glm")) {
      message0("Testing varHom ... ")
      FV.Model <- tryCatch(
        expr = update(F.Model, weights = NULL),
        warning = function(war) {
          message0("* Testing VarHom failed with the warning (see below): ")
          message0("\t", war, breakLine = FALSE)
          return(NULL)
        },
        error = function(err) {
          message0("* Testing VarHom failed with the error (see below): ")
          message0("\t", err, breakLine = FALSE)
          return(NULL)
        }
      )
      ########
      FV.Model <- intervalsCon(object = FV.Model, lvls = ci_levels)
      ######## following two lines are important
      FV.ML <- REML2ML(FV.Model)
      F.ML <- REML2ML(F.Model)
      if (!is.null(F.Model) &&
        !is.null(FV.Model) &&
        !is.null(FV.ML) &&
        !is.null(F.ML) &&
        # k is useless here
        AICc(FV.ML, k = log(n)) <= AICc(F.ML, k = log(n))) {
        # = is important
        F.Model <- FV.Model
        VarHomo <- FALSE
        message0("\tVarHom checked out ... ")
      }
    } else {
      optimise[2] <- FALSE
    }
    ###########
    if (optimise[3] && Batch_exist && !(mdl %in% "glm")) {
      message0("Testing Batch ... ")
      G.Model <- tryCatch(
        expr = do.call("gls",
          args =
            listFun(
              list = list(
                model     = formula(F.Model),
                data      = data,
                na.action = na.omit,
                method    = "REML",
                weights   = F.Model$call$weights,
                control   = gCont,
                ...
              ),
              FUN = "gls",
              debug = FALSE
            )
        ),
        warning = function(war) {
          message0("* Testing Batch failed with the warning (see below): ")
          message0("\t", war, breakLine = FALSE)
          return(NULL)
        },
        error = function(err) {
          message0("* Testing Batch failed with the error (see below): ")
          message0("\t", err, breakLine = FALSE)
          return(NULL)
        }
      )
      ########
      G.Model <- intervalsCon(object = G.Model, lvls = ci_levels)
      ######## following two lines are important
      G.ML <- REML2ML(G.Model)
      F.ML <- REML2ML(F.Model)
      if (!is.null(G.Model) &&
        !is.null(F.Model) &&
        !is.null(G.ML) &&
        !is.null(F.ML) &&
        # k is useless here
        AICc(G.ML, k = log(n)) <= AICc(F.ML, k = log(n))) {
        F.Model <- G.Model
        message0("\tBatch checked out ... ")
      }
    } else {
      optimise[3] <- FALSE
    }
    ###########
    if (optimise[4]) {
      SplitModels <- SplitEffect(
        finalformula = formula(F.Model),
        fullModelFormula = fixed,
        F.Model = F.Model,
        data = data,
        depVariable = allVars[1],
        ci_levels = ci_levels
      )
    }
    ###########
    if (optimise[5]) {
      message0("Estimating effect sizes ... ")
      EffectSizes <- c(suppressMessages(
        AllEffSizes(
          object = F.Model,
          depVariable = allVars[1],
          effOfInd = allVars[-1],
          data = data
        )
      ),
      "Combined effect sizes" = suppressMessages(if (!is.null(SplitModels)) {
        lapply(SplitModels, function(x) {
          percentageChangeCont(
            model = x,
            data = getData(x),
            variable = NULL,
            depVar = allVars[1],
            individual = FALSE,
            mainEffsOnlyWhenIndivi = x$MainEffect
          )
        })
      } else {
        NULL
      })
      )
      message0(
        "\tTotal effect sizes estimated: ",
        ifelse(
          !is.null(EffectSizes),
          length(EffectSizes),
          "Not possible because of errors!"
        )
      )
    }
    ###########
    if (optimise[6]) {
      message0("Quality tests in progress ... ")
      ResidualNormalityTest <- QuyalityTests(F.Model, levels = allVars[-1])
    }
  } else {
    message0("This process fully terminated. No success in recovering the initial model.")
    OutR$messages <- "This process fully terminated. No success in recovering the initial model."
    return(OutR)
  }
  message0("MM framework executed in ", round(difftime(Sys.time(), sta.time, units = "sec"), 2), " second(s).")
  ####
  OutR <- list(
    output = list(
      Final.Model = F.Model,
      Initial.Model = I.Model,
      "Effect sizes" = EffectSizes,
      ResidualNormalityTests = ResidualNormalityTest,
      NoVarStr.Model = FV.Model,
      NoBatch.Model = G.Model,
      SplitModels = SplitModels,
      Final.Model.Tag = class(F.Model)[1],
      Initial.Model.Tag = mdl,
      VarHomoIn = VarHomo,
      BatchIn = Batch_exist && is(F.Model, "lme"),
      SexIn = termInTheModel(
        model = formula(F.Model),
        term = "Sex",
        message = FALSE
      ),
      LifeStageIn = termInTheModel(
        model = formula(F.Model),
        term = "LifeStage",
        message = FALSE
      ),
      optimised = optimise
    ),
    input = list(
      OpenStatsList = object,
      fixed = initialFixed,
      random = random,
      data = data,
      depVariable = allVars[1],
      lower = lower,
      direction = direction,
      LifeStage = LifeStage,
      method = method,
      weight = weight,
      checks = checks,
      optimise = optimise,
      ci_level = ci_levels,
      others = ...
    ),
    extra = list(
      Cleanedformula = fixed,
      lowerCorrected = lowerCorrected
    )
  )
  class(OutR) <- "OpenStatsMM"
  return(OutR)
}
