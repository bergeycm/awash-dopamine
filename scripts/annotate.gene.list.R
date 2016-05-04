#!/usr/bin/env Rscript

annotate.gene.list = function(	gene.list, 
								uniprot.col.name="uniprot", 
								which.anno = c(	"pathway", "component", "go.mf", "go.bp", 
												"go.cc", "panther.proteinclass")) {

	# --------------------------------------------------------------------------------
	# --- Grab annotations of this type from PANTHER for these genes
	# --------------------------------------------------------------------------------

	anno.est = lapply(gene.list[[uniprot.col.name]], 
						function (x) {
							data.frame(gene=x, annotation=query.panther(x, which.anno))
						})
	anno.est = do.call(rbind, anno.est)

	# --------------------------------------------------------------------------------
	# --- Combine with original gene list, which may contain smoothed estimates
	# --------------------------------------------------------------------------------

	anno.gene.list = merge(	x=gene.list, y=anno.est, 
							by.x="uniprot", by.y="gene")

	return(anno.gene.list)
}
