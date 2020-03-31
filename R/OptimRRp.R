# Reference range OPT core
RRrunner = function(object              ,
										formula       = ~ category + Genotype + Sex + LifeStage,
										rep           = 1500                                   ,
										method        = NULL                                   ,
										RRrefLevel    = NULL                                   ,
										RRprop        = .05                                    ,
										ci_levels     = 0.95                                   ,
										fullComparisions = TRUE                                ,
										InterLevelComparisions = TRUE                          ,
										...)
{
	requireNamespace("rlist")
	sta.time    = Sys.time()
	allTerms    = all_vars0(x = formula)
	if (!method %in% c('RR')          ||
			is.null(allTerms)             ||
			is.null(object)               ||
			sum(allTerms %in% names(object@datasetPL)) < 2) {
		#####
		message0 (
			'Improper method (',
			method,
			') for the type of data, or the `formula/data` is not properly specified/left blank. \n\tFormula: ',
			printformula(formula)
		)
		return(NULL)
	}
	#message0('RR+ framework in progress ...')
	if (is.null(RRprop)     ||
			!is.numeric(RRprop) ||
			RRprop <= 0.5       ||
			RRprop >= 1.0) {
		message0('`RRprop` must be a value greater than 0.5 and less than 1')
		warnings('Improper value for "RRprop"')
		return(NULL)
	}
	RRpropTrans             = MakeRRQuantileFromTheValue(RRprop)
	cleanFormulaForOutput   = checkModelTermsInData(
		formula            = formula          ,
		data               = object@datasetPL ,
		responseIsTheFirst = TRUE
	)
	message0('Discritizing the continuous data into discrete levels. The quantile = ',
					 RRpropTrans)
	message0('Stp 1. Low versus Normal/High')
	RRobject_low = RRDiscretizedEngine(
		data     = object@datasetPL                   ,
		formula  = cleanFormulaForOutput              ,
		depVar   = allTerms[1]                        ,
		lower    = allTerms[2]                        ,
		refLevel = RRrefLevel                         ,
		labels   = c('Low', 'NormalHigh')             ,
		depVarPrefix = 'Low'                          ,
		right    = TRUE                               ,
		prob     = 1 - RRpropTrans
	)
	message0('Stp 2. Low/Normal versus High')
	RRobject_high = RRDiscretizedEngine(
		data     = object@datasetPL                   ,
		formula  = cleanFormulaForOutput              ,
		depVar   = allTerms[1]                        ,
		lower    = allTerms[2]                        ,
		refLevel = RRrefLevel                         ,
		labels   = c('LowNormal', 'High')             ,
		depVarPrefix = 'High'                         ,
		right    = FALSE                              ,
		prob     = RRpropTrans
	)
	###########################
	message0('Fisher exact test with '         ,
					 ifelse(rep > 0, rep, 'No')        ,
					 ' iteration(s) in progress ...')
	message0('Analysing Low vs NormalHigh ...')
	RRresult_low = lapply(RRobject_low, function(x) {
		r = suppressMessages(
			crunner(
				object   = x$newobject                       ,
				formula  = x$newFormula                      ,
				rep      = rep                               ,
				method   = 'RR'                              ,
				fullComparisions = fullComparisions          ,
				InterLevelComparisions = InterLevelComparisions ,
				noteToFinish     = 'in Low vs NormalHigh'    ,
				ci_levels        = ci_levels                 ,
				RRextraResults   = list(   
					depVariable      = x$depVariable           ,
					disdepVariable   = x$newDepVariable        ,
					RRpropTransformed= x$RRprop                ,
					RRLabels         = x$labels                ,
					RRprefix         = x$depVarPrefix          ,
					RRreferenceLevel = x$refLevel              ,
					RRempiricalQuantiles  = x$empiricalQuantiles
				)                                            ,
				trimWC = FALSE                               ,
				...
			)
		)
		return(r$output$SplitModels)
	})
	message0('Analysing LowNormal vs High ...')
	RRresult_high = lapply(RRobject_high, function(x) {
		r = suppressMessages(
			crunner(
				object   = x$newobject                       ,
				formula  = x$newFormula                      ,
				rep      = rep                               ,
				method   = 'RR'                              ,
				fullComparisions = fullComparisions          ,
				InterLevelComparisions = InterLevelComparisions ,
				noteToFinish     = 'in LowNormal vs High'    ,
				ci_levels        = ci_levels                 ,
				RRextraResults   = list(   
					depVariable      = x$depVariable           ,
					disdepVariable   = x$newDepVariable        ,
					RRpropTransformed= x$RRprop                ,
					RRLabels         = x$labels                ,
					RRprefix         = x$depVarPrefix          ,
					RRreferenceLevel = x$refLevel              ,
					RRempiricalQuantiles  = x$empiricalQuantiles
				)                                            ,
				trimWC             = FALSE                   ,  
				...
			)
		)
		return(r$output$SplitModels)
	})
	message0('RR framework executed in ', round(difftime(Sys.time() , sta.time, units = 'sec'), 2), ' second(s).')
	#####
	SpltResult = lapply(c(Low  = RRresult_low,
												High = RRresult_high), function(x) {
													x[[1]]
												})
	OutR = list(
		output = list(SplitModels = list.clean(SpltResult)) ,
		input  = list(
			OpenStatsList    = object                         ,
			data            = object@datasetPL                ,
			depVariable     = allTerms[1]                     ,
			rep             = rep                             ,
			method          = method                          ,
			formula         = formula                         ,
			prop            = RRprop                          ,
			ci_level        = ci_levels                       ,
			refLevel        = RRrefLevel                      ,
			full_comparisions = c(fullComparisions, InterLevelComparisions)
		),
		extra  = list(Cleanedformula           = cleanFormulaForOutput,
									TransformedRRprop        = RRpropTrans)
	)
	class(OutR) <- 'OpenStatsRR'
	return(OutR)
}
