sexbiasedgene = ddply(.data = T8, .variables = .(geneID, event_type), .fun = summarize, sexbiased = any(grepl("m|f", DevelopmentalSexBias)), .progress = 'text')
sexbiasedgene2 = recast(data = sexbiasedgene, formula = geneID~event_type)
