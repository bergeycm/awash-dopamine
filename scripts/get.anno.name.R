#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------------
# --- Get names of annotations (pathway names, cell components, etc.) given IDs
# ------------------------------------------------------------------------------------

get.anno.name = function (anno.id, all.anno.vals) {
	anno.name = strsplit(grep(paste0(anno.id, "$"), all.anno.vals, value=TRUE), 
							split="#")[[1]][1]
	return(anno.name)
}
