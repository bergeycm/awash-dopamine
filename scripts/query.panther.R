#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------------
# --- Query PANTHER data by Uniprot accession
# ------------------------------------------------------------------------------------

query.panther = function(acc, to.return=c("pathway", "component", "go.mf", "go.bp", 
							"go.cc", "panther.proteinclass")) {

	idx.uni = which(panther.classification$data$uniprot.id == acc)
	
	# Abort if not found in PANTHER
	if (length(idx.uni) == 0) {
		return(NA)
	}
	
	# Pathway or component
	if ((to.return == "pathway") || (to.return == "component")) {
	
		# Ensure length of panther.classification$panther.pathway[idx.uni] is 1
		if (length(panther.classification$panther.pathway[idx.uni]) > 1) {
			warning("Uniprot accession found more than once in PANTHER.")
		}
		
		query.result = panther.classification$panther.pathway[idx.uni][[1]]
		
		if (length(query.result) == 0) {
			return (NA)
		}
		
		# Pathways are those odd values that end with ">"
		if (to.return == "pathway") {
			result.clean = gsub("[>#]", "", grep(">", query.result, 
								value=TRUE))
		} else {
			result.clean = gsub("#", "", grep(">", query.result, 
								value=TRUE, invert=TRUE))
		}
	
	# GO terms and protein classes
	} else if (to.return %in% c("go.mf", "go.bp", "go.cc", "panther.proteinclass")) {

		# Ensure length of panther.classification$panther.pathway[idx.uni] is 1
		if (length(panther.classification[[to.return]][idx.uni]) > 1) {
			warning("Uniprot accession found more than once in PANTHER.")
		}
		
		query.result = panther.classification[[to.return]][idx.uni][[1]]

		if (length(query.result) == 0) {
			return (NA)
		}
		
		# Strip "GO:" or "#" from values
		result.clean = gsub("^#", "", gsub("^GO:", "", query.result))
		 
	}
	
	return(result.clean)
}
