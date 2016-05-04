#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------------
# --- Do enrichment test for one annotation value
# ------------------------------------------------------------------------------------

enrich.test.annotation = function (annotated.gene.list, anno.val, ref.list) {

	# Do Wilcox test
	
	matches.nodup = get.anno.matches(annotated.gene.list, anno.val)
	n.matches = nrow(matches.nodup)

	# Abort if no matches
	if (n.matches == 0) {
		return(c(0, NA, NA, NA))
	}

	wilcox = wilcox.test(	x=matches.nodup$estimate, 
							y=ref.list$estimate)

	if (mean(matches.nodup$estimate) > mean(ref.list$estimate)) {
		direction = "+"
	} else {
		direction = "-"
	}
		
	return(c(n.matches, wilcox$statistic, direction, as.numeric(wilcox$p.value)))

}
