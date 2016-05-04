#!/usr/bin/env Rscript

# --------------------------------------------------------------------------------
# --- Do enrichment test for each annotation value
# --- (each pathway, each cell component, etc.)
# --------------------------------------------------------------------------------

enrich.test.all = function (annotated.gene.list, all.annotations) {

	# Get list of all genes to use as comparison in enrichment tests
	all = annotated.gene.list[!duplicated(annotated.gene.list$gene.cap),]

	# Bring in function to do enrichment test for single annotation value
	source("scripts/enrich.test.annotation.R")

	anno.enrich.results = as.data.frame(t(sapply(	all.annotations, 
													enrich.test.annotation,
													annotated.gene.list = annotated.gene.list,
													ref.list=all)), 
										stringsAsFactors=FALSE)
	names(anno.enrich.results) = c("N", "W", "overUnder", "p")

	return(anno.enrich.results)
}
