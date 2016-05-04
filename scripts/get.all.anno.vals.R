#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------------
# --- Get all possible values for a given annotation type
# ------------------------------------------------------------------------------------

get.all.anno.vals = function(this.anno, panther) {

	if (this.anno == "pathway") {
		anno.split = unique(unlist(strsplit(panther$panther.pathway, split="[>;]")))
	} else {
		anno.split = unique(unlist(strsplit(panther[[this.anno]], split=";")))
	}
	
	return(anno.split)
}
