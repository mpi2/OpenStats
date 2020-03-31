OpenStatsReportCont =	function(object,
															 debug = FALSE)
{
	if (!is.null(object$messages))
		return (NULL)
	#####################################################################
	Labels         = OpenStatsListLevels(object = object)
	Fmodel         = object$output$Final.Model
	frm            = formula(Fmodel)
	fim            = object$input$fixed
	depVariable    = all_vars0(frm)[1]
	equation       = ifelse(
		Labels$Weight %in% all_vars0(frm),
		paste0('including '    , Labels$Weight),
		paste0('not including ', Labels$Weight)
	)
	formula        = printformula(frm)
	#modelContrast  = modelContrasts(formula = frm,data = object$input$data)
	framework      = switch(
		#object$output$Final.Model.Tag 
		class(object$output$Final.Model),
		lme  = "Linear Mixed Model framework"                          ,
		gls  = "Linear Model Using Generalized Least Squares framework",
		glm  = "Generalized Linear Model framework"
	)
	#fittingMethod    = toupper(object$output$Final.Model.Tag)
	fittingMethod    = toupper(class(object$output$Final.Model))
	#####################################################################
	x                = object$input$OpenStatsList@datasetPL
	columnOfInterest = x[, c(depVariable)]
	#####################################################################
	variability      = list('Value' = length(unique(columnOfInterest)) / max(length(columnOfInterest), 1), 
													'Type'  = 'Total unique response divided by total number of response')
	#####################################################################
	DSsize            = SummaryStats(
		x = object$input$data         ,
		formula = checkModelTermsInData(
			formula = object$input$fixed,
			data = x,
			responseIsTheFirst = TRUE
		),
		#label = 'Summary statistics',
		lower  = TRUE,
		drop   = TRUE,
		sep    = '_'
	)
	MultiBatch = ifelse(multiBatch(x),
											'Dataset contains multiple batches',
											'Dataset contains single batch')
	addInfo           = list(
		Data = list(
			'Data signature'         = dataSignature(formula = object$input$fixed,
																							 data    = object$input$data),
			'Variability'            = variability                               ,
			'Summary statistics'     = DSsize
		),
		# 'Formula'                = list(
		# 	input   = printformula(object$input$fixed),
		# 	final   =	printformula(formula)
		# ),
		Analysis = list(
			'Model setting'          = extractLmeTerms(object)         ,
			'Is model optimised'     = optimM(object$output$optimised) ,
			'Multibatch in analysis' = MultiBatch                      ,
			'Gender included in analysis' = GenderIncludedInAnalysis(x),
			'Further models' = if (!is.null(object$output$SplitModels)) {
				lapply(object$output$SplitModels, function(v) {
					if (class(v) %in% c('lme', 'gls', 'glm')) {
						r          = as.list(unmatrix0(summary1(v)$tTable))
						r$Model    = printformula(v$SplitFormula)
						r$Method   = pasteComma(class(v),truncate = FALSE)
					} else{
						r = v
					}
					return(r)
				})
			} else{
				NULL
			},
			'Effect sizes'                   = object$output$'Effect sizes',
			'Other residual normality tests' = object$output$ResidualNormalityTests
		)
	)
	#####################################################################
	pcS = object$output$'Effect sizes'$'Combined effect sizes.Genotype_Sex'$'Percentage change'
	pcO = object$output$'Effect sizes'$Genotype$'Percentage change'
	percentageChanges = if (!is.null(pcS)) {
		pcS
	} else{
		pcO
	}
	SexDymFinalModel = TermInFormulaReturn(
		active = TRUE,
		formula = frm,
		term = CombineLevels(Labels$Genotype$Genotype, Labels$Sex$Sex, debug = debug),
		return = TRUE,
		not = FALSE,
		debug = debug
	)
	SexDymInputModel = TermInFormulaReturn(
		active = TRUE,
		formula = fim,
		term = CombineLevels(Labels$Genotype$Genotype, Labels$Sex$Sex, debug = debug),
		return = TRUE,
		not = FALSE,
		debug = debug
	)
	
	#####################################################################
	OpenStatsReportMM0      = list(
		'Applied method'                       = 	paste0(framework, ", ", fittingMethod, ', ', format(equation)),
		'Dependent variable'                   =	depVariable              ,
		'Batch included'                       =	object$output$BatchIn    ,
		'Batch p-value'                        =  NULL                     ,
		'Residual variances homogeneity'       =	object$output$VarHomoIn  ,
		'Residual variances homogeneity p-value' =  NULL                     ,
		#####################################################################
		'Genotype contribution' =	list(
			Overall = TermInFormulaReturn(
				active = TRUE,
				formula = frm,
				term = CombineLevels(Labels$Genotype$Genotype,
														 Labels$Sex$Sex,
														 debug = debug),
				return = NULL,
				not = modelSummaryPvalueExtract(
					x = Fmodel,
					variable = Labels$Genotype$Genotype,
					anova = TRUE,
					debug = debug
				),
				debug = debug
			),
			'Sex FvKO p-value'   =	TermInFormulaReturn(
				active = TRUE,
				formula = frm,
				term = CombineLevels(
					Labels$Sex$Sex,
					Labels$Genotype$Genotype,
					debug = debug
				),
				not = NULL,
				return = modelSummaryPvalueExtract(
					x = object$output$SplitModels$Genotype_Sex,
					# SexFemale:Genotypeexperimental
					variable = CombineLevels(
						paste0('Sex', Labels$Sex$Female),
						Labels$Genotype$Levels,
						debug = debug
					),
					anova = FALSE,
					debug = debug
				),
				debug = debug
			),
			'Sex MvKO p-value'  =	TermInFormulaReturn(
				active = TRUE,
				formula = frm,
				term = CombineLevels(
					Labels$Sex$Sex,
					Labels$Genotype$Genotype,
					debug = debug
				),
				not = NULL,
				return = modelSummaryPvalueExtract(
					x = object$output$SplitModels$Genotype_Sex,
					variable = CombineLevels(
						paste0('Sex', Labels$Sex$Male),
						Labels$Genotype$Levels,
						debug = debug
					),
					anova = FALSE,
					debug = debug
				),
				debug = debug
			),
			'Sexual dimorphism detected' = list(
				'Criteria' = SexDymFinalModel              ,
				'Note'     = paste0(
					'Genotype-Sex interaction '              ,
					ifelse(SexDymInputModel	,	'is',	'is not'),
					' part of the input '                    ,
					ifelse(SexDymFinalModel,	'(it is part of the final)',	'(it is not part of the final)'),
					' model. '
				)
			)
		),
		'Genotype estimate' =
			modelSummaryPvalueExtract(
				x = Fmodel,
				variable = unlist(Labels$Genotype$Levels),
				anova = FALSE,
				what = 'Value',
				debug = debug,
				ci_display = TRUE
			),
		'Genotype standard error'  =
			modelSummaryPvalueExtract(
				x = Fmodel,
				variable = unlist(Labels$Genotype$Levels),
				anova = FALSE,
				what = 'Std.Error',
				debug = debug
			),
		'Genotype p-value'           =
			modelSummaryPvalueExtract(
				x = Fmodel,
				variable = Labels$Genotype$Genotype,
				anova = TRUE,
				debug = debug
			),
		'Genotype percentage change'           =	percentageChanges,
		'Genotype effect size'                 = object$output$'Effect sizes'[[Labels$Genotype$Genotype]],
		#####################################################################
		'Sex estimate'                         =	modelSummaryPvalueExtract(
			x = Fmodel,
			variable = unlist(Labels$Sex$Levels),
			anova = FALSE,
			what = 'Value',
			debug = debug,
			ci_display = TRUE
		),
		'Sex standard error'                   = modelSummaryPvalueExtract(
			x = Fmodel,
			variable = unlist(Labels$Sex$Levels),
			anova = FALSE,
			what = 'Std.Error',
			debug = debug
		),
		'Sex p-value'                            =	modelSummaryPvalueExtract(
			x = Fmodel,
			variable = Labels$Sex$Sex,
			anova = TRUE,
			debug = debug
		),
		'Sex effect size'                         =	object$output$'Effect sizes'[[Labels$Sex$Sex]],
		#####################################################################
		'LifeStage estimate'                      =	modelSummaryPvalueExtract(
			x = Fmodel,
			variable = unlist(Labels$LifeStage$Levels),
			anova = FALSE,
			what = 'Value',
			debug = debug,
			ci_display = TRUE
		),
		'LifeStage standard error'                =	modelSummaryPvalueExtract(
			x = Fmodel,
			variable = unlist(Labels$LifeStage$Levels),
			anova = FALSE,
			what = 'Std.Error',
			debug = debug
		),
		'LifeStage p-value'                         =	modelSummaryPvalueExtract(
			x = Fmodel,
			variable = Labels$LifeStage$LifeStage,
			anova = TRUE,
			debug = debug
		),
		'LifeStage effect size'                = object$output$'Effect sizes'[[Labels$LifeStage$LifeStage]] ,
		#####################################################################
		'Weight estimate'                      =	modelSummaryPvalueExtract(
			x = Fmodel,
			variable = Labels$Weight,
			anova = FALSE,
			what = 'Value',
			debug = debug,
			ci_display = TRUE
		),
		'Weight standard error'                =	modelSummaryPvalueExtract(
			x = Fmodel,
			variable = Labels$Weight,
			anova = FALSE,
			what = 'Std.Error',
			debug = debug
		),
		'Weight p-value'                         =	modelSummaryPvalueExtract(
			x = Fmodel,
			variable = Labels$Weight,
			anova = TRUE,
			debug = debug
		),
		'Weight effect size'                   =  object$output$'Effect sizes'[[Labels$Weight]],
		#####################################################################
		'Gp1 genotype'                         =	Labels$Genotype$Control		,
		'Gp1 Residuals normality test'         =	object$output$ResidualNormalityTests$Genotype[Labels$Genotype$Control][[1]],
		'Gp2 genotype'                         =	Labels$Genotype$Mutant				,
		'Gp2 Residuals normality test'         =	object$output$ResidualNormalityTests$Genotype[Labels$Genotype$Mutant][[1]],
		#####################################################################
		'Blups test'                           =  NULL,
		'Rotated residuals normality test'     =  NULL,
		#####################################################################
		'Intercept estimate'                   =	modelSummaryPvalueExtract(
			x = Fmodel,
			variable = '(Intercept)',
			anova = FALSE,
			what = 'Value',
			debug = debug,
			ci_display = TRUE
		),
		'Intercept standard error'             =	modelSummaryPvalueExtract(
			x = Fmodel,
			variable = '(Intercept)',
			anova = FALSE,
			what = 'Std.Error',
			debug = debug
		),
		'Intercept p-value'                      =	modelSummaryPvalueExtract(
			x = Fmodel,
			variable = '(Intercept)',
			anova = TRUE,
			debug = debug
		),
		#####################################################################
		'Interactions included'              =	list(
			'Genotype Sex'  =  TermInFormulaReturn(
				formula = frm,
				term = CombineLevels(Labels$Genotype$Genotype, Labels$Sex$Sex, debug = debug),
				return =
					!is.null(
						modelSummaryPvalueExtract(
							x = Fmodel,
							variable = CombineLevels(Labels$Genotype$Genotype, Labels$Sex$Sex, debug = debug),
							anova = TRUE,
							debug = debug
						)
					),
				not = NULL
			),
			'Genotype LifeStage'  =  TermInFormulaReturn(
				formula = frm,
				term = CombineLevels(Labels$Genotype$Genotype, Labels$LifeStage$LifeStage, debug = debug),
				return = !is.null(
					modelSummaryPvalueExtract(
						x = Fmodel,
						variable = CombineLevels(Labels$Genotype$Genotype, Labels$LifeStage$LifeStage, debug = debug),
						anova = TRUE,
						debug = debug
					)
				),
				not = NULL
			)
			,
			'Sex LifeStage'  =  TermInFormulaReturn(
				formula = frm,
				term = CombineLevels(Labels$Sex$Sex, Labels$LifeStage$LifeStage, debug = debug),
				return = !is.null(
					modelSummaryPvalueExtract(
						x = Fmodel,
						variable = CombineLevels(Labels$Sex$Sex, Labels$LifeStage$LifeStage, debug = debug),
						anova = TRUE,
						debug = debug
					)
				),
				not = NULL
			),
			'Genotype Sex LifeStage'  =  TermInFormulaReturn(
				formula = frm,
				term = CombineLevels(
					Labels$Genotype$Genotype  ,
					Labels$Sex$Sex            ,
					Labels$LifeStage$LifeStage,
					len   = 3                 ,
					debug = debug
				),
				return = !is.null(
					modelSummaryPvalueExtract(
						x = Fmodel,
						variable = CombineLevels(
							Labels$Genotype$Genotype  ,
							Labels$Sex$Sex            ,
							Labels$LifeStage$LifeStage,
							len   = 3                 ,
							debug = debug
						),
						anova = TRUE,
						debug = debug
					)
				),
				not = NULL
			)
		),
		#####################################################################
		################ interaction
		'Interactions p-value'      =	list(
			'Genotype Sex'  =
				modelSummaryPvalueExtract(
					x = Fmodel,
					variable =  CombineLevels(Labels$Genotype$Genotype, Labels$Sex$Sex,debug = debug),
					anova = TRUE,
					debug = debug
				),
			'Genotype LifeStage'  =
				modelSummaryPvalueExtract(
					x = Fmodel,
					variable =  CombineLevels(Labels$Genotype$Genotype, Labels$LifeStage$LifeStage,debug = debug),
					anova = TRUE,
					debug = debug
				),
			'Sex LifeStage'  =
				modelSummaryPvalueExtract(
					x = Fmodel,
					variable =  CombineLevels(Labels$Sex$Sex, Labels$LifeStage$LifeStage,debug = debug),
					anova = TRUE,
					debug = debug
				),
			'Genotype Sex LifeStage'  =
				modelSummaryPvalueExtract(
					x = Fmodel,
					variable =  CombineLevels(
						Labels$Genotype$Genotype   ,
						Labels$Sex$Sex             ,
						Labels$LifeStage$LifeStage ,
						debug = debug              ,
						len = 3
					),
					anova = TRUE,
					debug = debug
				)
		),
		#####################################################################
		################ Sex interactions
		'Sex FvKO estimate'                    =
			modelSummaryPvalueExtract(
				x = object$output$SplitModels$Genotype_Sex,
				variable = CombineLevels(
					paste0('Sex', Labels$Sex$Female),
					Labels$Genotype$Levels,
					debug = debug
				),
				anova = FALSE,
				what = 'Value',
				debug = debug,
				ci_display = TRUE
			),
		'Sex FvKO standard error'              =	modelSummaryPvalueExtract(
			x = object$output$SplitModels$Genotype_Sex,
			variable = CombineLevels(
				paste0('Sex', Labels$Sex$Female),
				Labels$Genotype$Levels,
				debug = debug
			),
			anova = FALSE,
			what = 'Std.Error',
			debug = debug
		),
		'Sex FvKO p-value'                       = modelSummaryPvalueExtract(
			x = object$output$SplitModels$Genotype_Sex,
			variable = CombineLevels(
				paste0('Sex', Labels$Sex$Female),
				Labels$Genotype$Levels,
				debug = debug
			),
			anova = FALSE,
			debug = debug
		),
		'Sex FvKO effect size'                 = object$output$'Effect sizes'[[paste(Labels$Genotype$Genotype, Labels$Sex$Female, sep = '_')]],
		#####################################################################
		'Sex MvKO estimate'                    =	modelSummaryPvalueExtract(
			x = object$output$SplitModels$Genotype_Sex,
			variable = CombineLevels(
				paste0('Sex', Labels$Sex$Male),
				Labels$Genotype$Levels,
				debug = debug
			),
			anova = FALSE,
			what = 'Value',
			debug = debug,
			ci_display = TRUE
		),
		'Sex MvKO standard error'              =	 modelSummaryPvalueExtract(
			x = object$output$SplitModels$Genotype_Sex,
			variable = CombineLevels(
				paste0('Sex', Labels$Sex$Male),
				Labels$Genotype$Levels,
				debug = debug
			),
			anova = FALSE,
			what = 'Std.Error',
			debug = debug
		),
		'Sex MvKO p-value'                       =	modelSummaryPvalueExtract(
			x = object$output$SplitModels$Genotype_Sex,
			variable = CombineLevels(
				paste0('Sex', Labels$Sex$Male),
				Labels$Genotype$Levels,
				debug = debug
			),
			anova = FALSE,
			debug = debug
		),
		'Sex MvKO effect size'                 = object$output$'Effect sizes'[[paste(Labels$Genotype$Genotype, Labels$Sex$Male, sep = '_')]],
		#####################################################################
		################ LifeStage interaction
		'LifeStage EvKO estimate'                    =	modelSummaryPvalueExtract(
			x = object$output$SplitModels$Genotype_LifeStage,
			variable = CombineLevels(
				paste0(Labels$LifeStage$LifeStage, Labels$LifeStage$Early),
				Labels$Genotype$Levels,
				debug = debug
			),
			anova = FALSE,
			what = 'Value',
			debug = debug,
			ci_display = TRUE
		),
		'LifeStage EvKO standard error'              =	modelSummaryPvalueExtract(
			x = object$output$SplitModels$Genotype_LifeStage,
			variable =  CombineLevels(
				paste0(Labels$LifeStage$LifeStage, Labels$LifeStage$Early),
				Labels$Genotype$Levels,
				debug = debug
			),
			anova = FALSE,
			what = 'Std.Error',
			debug = debug
		),
		'LifeStage EvKO p-value'                       =	modelSummaryPvalueExtract(
			x = object$output$SplitModels$Genotype_LifeStage,
			variable =  CombineLevels(
				paste0(Labels$LifeStage$LifeStage, Labels$LifeStage$Early),
				Labels$Genotype$Levels,
				debug = debug
			),
			anova = FALSE,
			debug = debug
		),
		'LifeStage EvKO effect size'                 = object$output$'Effect sizes'$Genotype_Early ,
		#####################################################################
		'LifeStage LvKO estimate'                    =	modelSummaryPvalueExtract(
			x = object$output$SplitModels$Genotype_LifeStage,
			variable =  CombineLevels(
				paste0(Labels$LifeStage$LifeStage, Labels$LifeStage$Late),
				Labels$Genotype$Levels,
				debug = debug
			),
			anova = FALSE,
			what = 'Value',
			debug = debug,
			ci_display = TRUE
		),
		'LifeStage LvKO standard error'              =	modelSummaryPvalueExtract(
			x = object$output$SplitModels$Genotype_LifeStage,
			variable = CombineLevels(
				paste0(Labels$LifeStage$LifeStage, Labels$LifeStage$Late),
				Labels$Genotype$Levels,
				debug = debug
			),
			anova = FALSE,
			what = 'Std.Error',
			debug = debug
		),
		'LifeStage LvKO p-value'                       =	modelSummaryPvalueExtract(
			x = object$output$SplitModels$Genotype_LifeStage,
			variable = CombineLevels(
				paste0(Labels$LifeStage$LifeStage, Labels$LifeStage$Late),
				Labels$Genotype$Levels,
				debug = debug
			),
			anova = FALSE,
			debug = debug
		),
		'LifeStage LvKO effect size'                 = object$output$'Effect sizes'$Genotype_Late ,
		#####################################################################
		################ Sex LifeStage Genotype interactions
		# 1.
		'LifeStageSexGenotype FvEvKO estimate'           =	modelSummaryPvalueExtract(
			x = object$output$SplitModels$Genotype_Sex.LifeStage,
			variable = 		CombineLevels(
				paste0(Labels$Sex$Sex, Labels$Sex$Female),
				paste0(Labels$LifeStage$LifeStage, Labels$LifeStage$Early),
				Labels$Genotype$Levels,
				debug = debug
			),
			anova = FALSE,
			what = 'Value',
			debug = debug,
			ci_display = TRUE
		),
		'LifeStageSexGenotype FvEvKO standard error'           =	modelSummaryPvalueExtract(
			x = object$output$SplitModels$Genotype_Sex.LifeStage,
			variable = 		CombineLevels(
				paste0(Labels$Sex$Sex, Labels$Sex$Female),
				paste0(Labels$LifeStage$LifeStage, Labels$LifeStage$Early),
				Labels$Genotype$Levels,
				debug = debug
			),
			anova = FALSE,
			what = 'Std.Error',
			debug = debug
		),
		'LifeStageSexGenotype FvEvKO p-value'           =	modelSummaryPvalueExtract(
			x = object$output$SplitModels$Genotype_Sex.LifeStage,
			variable = 		CombineLevels(
				paste0(Labels$Sex$Sex, Labels$Sex$Female),
				paste0(Labels$LifeStage$LifeStage, Labels$LifeStage$Early),
				Labels$Genotype$Levels,
				debug = debug
			),
			anova = FALSE,
			debug = debug
		),
		'LifeStageSexGenotype FvEvKO effect size'        = object$output$'Effect sizes'$'Genotype_Female Early' ,
		# 2.
		'LifeStageSexGenotype MvEvKO estimate'           =	modelSummaryPvalueExtract(
			x = object$output$SplitModels$Genotype_Sex.LifeStage,
			variable = 		CombineLevels(
				paste0(Labels$Sex$Sex, Labels$Sex$Male),
				paste0(Labels$LifeStage$LifeStage, Labels$LifeStage$Early),
				Labels$Genotype$Levels,
				debug = debug
			),
			anova = FALSE,
			what = 'Value',
			debug = debug,
			ci_display = TRUE
		),
		'LifeStageSexGenotype MvEvKO standard error'           =	modelSummaryPvalueExtract(
			x = object$output$SplitModels$Genotype_Sex.LifeStage,
			variable = CombineLevels(
				paste0(Labels$Sex$Sex, Labels$Sex$Male),
				paste0(Labels$LifeStage$LifeStage, Labels$LifeStage$Early),
				Labels$Genotype$Levels,
				debug = debug
			),
			anova = FALSE,
			what = 'Std.Error',
			debug = debug
		),
		'LifeStageSexGenotype MvEvKO p-value'           =	modelSummaryPvalueExtract(
			x = object$output$SplitModels$Genotype_Sex.LifeStage,
			variable = CombineLevels(
				paste0(Labels$Sex$Sex, Labels$Sex$Male),
				paste0(Labels$LifeStage$LifeStage, Labels$LifeStage$Early),
				Labels$Genotype$Levels,
				debug = debug
			),
			anova = FALSE,
			debug = debug
		),
		'LifeStageSexGenotype MvEvKO effect size'        = object$output$'Effect sizes'$'Genotype_Male Early' ,
		# 3.
		'LifeStageSexGenotype FvLvKO estimate'           =	modelSummaryPvalueExtract(
			x = object$output$SplitModels$Genotype_Sex.LifeStage,
			variable = CombineLevels(
				paste0(Labels$Sex$Sex, Labels$Sex$Female),
				paste0(Labels$LifeStage$LifeStage, Labels$LifeStage$Late),
				Labels$Genotype$Levels,
				debug = debug
			),
			anova = FALSE,
			what = 'Value',
			debug = debug,
			ci_display = TRUE
		),
		'LifeStageSexGenotype FvLvKO standard error'           =	modelSummaryPvalueExtract(
			x = object$output$SplitModels$Genotype_Sex.LifeStage,
			variable =  CombineLevels(
				paste0(Labels$Sex$Sex, Labels$Sex$Female),
				paste0(Labels$LifeStage$LifeStage, Labels$LifeStage$Late),
				Labels$Genotype$Levels,
				debug = debug
			),
			anova = FALSE,
			what = 'Std.Error',
			debug = debug
		),
		'LifeStageSexGenotype FvLvKO p-value'           =	modelSummaryPvalueExtract(
			x = object$output$SplitModels$Genotype_Sex.LifeStage,
			variable =  CombineLevels(
				paste0(Labels$Sex$Sex, Labels$Sex$Female),
				paste0(Labels$LifeStage$LifeStage, Labels$LifeStage$Late),
				Labels$Genotype$Levels,
				debug = debug
			),
			anova = FALSE,
			debug = debug
		),
		'LifeStageSexGenotype FvLvKO effect size'        = object$output$'Effect sizes'$'Genotype_Female Late' ,
		#4.
		'LifeStageSexGenotype MvLvKO estimate'           =	modelSummaryPvalueExtract(
			x = object$output$SplitModels$Genotype_Sex.LifeStage,
			variable =  CombineLevels(
				paste0(Labels$Sex$Sex, Labels$Sex$Male),
				paste0(Labels$LifeStage$LifeStage, Labels$LifeStage$Late),
				Labels$Genotype$Levels,
				debug = debug
			),
			anova = FALSE,
			what = 'Value',
			debug = debug,
			ci_display = TRUE
		),
		'LifeStageSexGenotype MvLvKO standard error'           =	modelSummaryPvalueExtract(
			x = object$output$SplitModels$Genotype_Sex.LifeStage,
			variable = CombineLevels(
				paste0(Labels$Sex$Sex, Labels$Sex$Male),
				paste0(Labels$LifeStage$LifeStage, Labels$LifeStage$Late),
				Labels$Genotype$Levels,
				debug = debug
			),
			anova = FALSE,
			what = 'Std.Error',
			debug = debug
		),
		'LifeStageSexGenotype MvLvKO p-value'           =	modelSummaryPvalueExtract(
			x = object$output$SplitModels$Genotype_Sex.LifeStage,
			variable = CombineLevels(
				paste0(Labels$Sex$Sex, Labels$Sex$Male),
				paste0(Labels$LifeStage$LifeStage, Labels$LifeStage$Late),
				Labels$Genotype$Levels,
				debug = debug
			),
			anova = FALSE,
			debug = debug
		),
		'LifeStageSexGenotype MvLvKO effect size'        = object$output$'Effect sizes'$'Genotype_Male Late' ,
		################
		'Classification tag'                   =	NULL,
		'Transformation'                       =	NULL,
		'Additional information'               =	addInfo
	)
	
	return(OpenStatsReportMM0)
}
