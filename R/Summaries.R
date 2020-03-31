summary.NULL = function(object, ...) {
	message0('No summary available for a NULL object')
}
summary.OpenStatsRR = function(object, format = 'rst', ...) {
	if (!is.null(object$messages) || is.null(object)) {
		message0('Due to error(s), no summary available')
		message0(object$messages)
		stop()
	}
	summaryCore(object, procedure = 'RR', format = format, ...)
}

summary.OpenStatsFE = function(object, format = 'rst', ...) {
	if (!is.null(object$messages) || is.null(object)) {
		message0('Due to error(s), no summary available')
		message0(object$messages)
		stop()
	}
	summaryCore(object, procedure = 'FE', format = format, ...)
}

summary.OpenStatsMM = function(object, format = 'rst', ...) {
	if (!is.null(object$messages) || is.null(object)) {
		message0('Due to error(s), no summary available')
		message0(object$messages)
		stop()
	}
	summaryCore(object, procedure = 'MM', format = format, ...)
}


summaryCore = function(x,
											 procedure = 'MM',
											 format = 'rst',
											 ...) {
	message0('Working on the summary table ...')
	requireNamespace("knitr")
	vo = OpenStatsReport(object           = x    ,
													JSON             = FALSE,
													ReportNullSchema = FALSE,
													RemoveNullKeys   = FALSE )
	pasteComma2 = function(...) {
		requireNamespace("rlist")
		inp = as.list(...)
		inp = list.clean(inp)
		if (is.null(inp)           ||
				length (inp) < 1)
			r = ' - '
		else
			r = pasteCommaJustForSummary(
				...,
				replaceNull   = TRUE  ,
				truncate      = FALSE ,
				replaceNullby = ' - '
			)
		return(r)
	}
	IfoverallTableExistsThenReturnOne = function(x, overallTableName  = 'Complete table') {
		if (overallTableName %in% names(x)) {
			#message0('`',overallTableName, '` found in the object')
			x = x[[overallTableName]]
		}
		return(x)
	}
	out = list(
		'Applied framework'              = vo$`Applied method`,
		'Final model'                    = if (procedure %in% 'MM') {
			if (!is.null(x$output$Final.Model))
				formula(x$output$Final.Model)
			else
				NULL
		} else if (procedure %in% 'FE') {
			RightFormula2LeftFormula(x$extra$Cleanedformula)
		}	else if (procedure %in% 'RR') {
			RightFormula2LeftFormula(x$extra$Cleanedformula)
		} else{
			NULL
		},
		'............................'   = '............................',
		'Tested Gene'                    = vo$`Gp2 genotype`,
		'Reference Gene'                 = vo$`Gp1 genotype`,
		'............................'   = '............................',
		'Sexual dimorphism detected?'    = pasteComma(vo$`Genotype contribution`$`Sexual dimorphism detected`),
		'............................'   = ifelse(
			procedure == 'RR',
			'* Separate p-values for (Low vs NormalHigh), (LowNormal vs High) and details ',
			'............................'
		),
		'Genotype contribution overall'  = pasteComma2(IfoverallTableExistsThenReturnOne(vo$`Genotype p-value`)),
		'Genotype contribution Females'  = pasteComma2(IfoverallTableExistsThenReturnOne(vo$`Sex FvKO p-value`)),
		'Genotype contribution Males'    = pasteComma2(IfoverallTableExistsThenReturnOne(vo$`Sex MvKO p-value`)),
		'............................'   = '............................'                                      ,
		'LifeStage contribution'         = pasteComma2(IfoverallTableExistsThenReturnOne(vo$`LifeStage p-value`)),
		'Genotype contribution Early'    = pasteComma2(IfoverallTableExistsThenReturnOne(vo$`LifeStage EvKO p-value`)),
		'Genotype contribution Late'     = pasteComma2(IfoverallTableExistsThenReturnOne(vo$`LifeStage LvKO p-value`)),
		'............................'   = '............................'                                      ,
		'Sex contribution'               = pasteComma2(IfoverallTableExistsThenReturnOne(vo$`Sex p-value`))     ,
		'Body weight contribution'       = pasteComma2(IfoverallTableExistsThenReturnOne(vo$`Weight p-value`))
	)
	outT = prepareSummaryOutput(out)
	print(kable(
		outT,
		format = format,
		col.names = c('Statistic', 'Value'),
		...
	))
	outTRemoved = outT[!apply(outT, 1, function(x) {
		any(grepl(
			pattern = '..........',
			x = x,
			fixed = TRUE
		))
	}), ]
	return(invisible(outTRemoved))
}

prepareSummaryOutput = function(out, nullMessage = 'Not applicable') {
	outT = as.matrix(out)
	outT = as.matrix(cbind(rownames(outT), outT))
	rownames(outT) = NULL
	return(outT)
}

pasteCommaJustForSummary = function(...,
											replaceNull   = TRUE   ,
											truncate      = TRUE   ,
											width         = 100    ,
											trailingSpace = TRUE   ,
											replaceNullby = 'NULL') {
	sep = ifelse(trailingSpace, ', ', ',')
	if (replaceNull)
		r = paste(
			replaceNull(as.list(...), replaceBy = replaceNullby),
			sep      = sep,
			collapse = sep
		)
	else
		r = paste(...,
							sep      = sep,
							collapse = sep)
	
	if (truncate)
		r = truncate_text(r, width)
	
	return(r)
}
