# --------------------------------------------------------------------------------
# --- Write this info (sans NA values) to file
# --------------------------------------------------------------------------------

write.overrep.test.results = function (	anno.overrep.results, est.name, 
										this.anno, anno.split) {

	anno.overrep.results.noNA = anno.overrep.results[
			sapply(anno.overrep.results, function(x){ sum(is.na(x[[1]])) == 0  })
	]
	
	anno.overrep.results.noNA = do.call(rbind, anno.overrep.results.noNA)
	anno.overrep.results.noNA$over.under = 
			unlist(anno.overrep.results.noNA$over.under)
	
	# Add columns with estimate and annotation value (e.g. fst and pathway)
	anno.overrep.results.noNA = cbind(est.name, this.anno, anno.overrep.results.noNA)
	
	# Add column with annotation name
	anno.overrep.results.noNA$name = sapply(anno.overrep.results.noNA$category,
											get.anno.name,
											all.anno.vals = anno.split)
	
	# Write table
	write.table(anno.overrep.results.noNA, 
				paste0(out.dir, est.name, ".overrep.", this.anno, ".txt"),
				row.names=FALSE)
	
}