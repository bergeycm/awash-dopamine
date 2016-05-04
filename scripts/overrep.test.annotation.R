#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------------
# --- Do overrepresentation test for one annotation value
# ------------------------------------------------------------------------------------

overrep.test.annotation = function (anno.val, est.name, annotated.gene.list, extremes) {

	write(paste0("Processing ", est.name, " and ", 
					this.anno, " - ", anno.val, "..."), stderr())
	
	matches.nodup = get.anno.matches(annotated.gene.list, anno.val)
	n.matches = nrow(matches.nodup)

	# Abort if no matches
	if (n.matches == 0) {
		return(list(rep(NA, 7)))
	}
	
	# Get list of all genes to use as reference
	all = annotated.gene.list[!duplicated(annotated.gene.list$gene.cap),]

	# Number of genes in reference list
	N = length(all$gene.cap)
	# Number of genes that map to category
	nC = n.matches
	# Based on this, percent of ref list genes that are in category
	pC = nC / N

	p.cutoffs = seq(from = 0.01, to = 0.10, by=0.005)

	overrep.results.list = list()
	
	for (dir in c("hi", "lo")) {
	
		overrep.results = as.data.frame(do.call(rbind, 
			lapply(p.cutoffs, function(x) {
				overrep.test.annotation.pval(x, high.or.low.p=dir, 
												extremes=extremes,
												all=all,
												matches.nodup=matches.nodup,
												pC=pC, nC=nC)
			})
		))
	
		overrep.results = cbind(anno.val, overrep.results, dir)
	
		names(overrep.results) = c("category", "p.val.cutoff", "cat.obs", 
									"outlier.cat.obs", "outlier.cat.exp", 
									"over.under", "fold.enrich",
									"overrep.p", "hi.or.lo.p.cutoff")

		# Unlist numeric columns. Hacky.
		num.cols = c(2:5, 7:8)
		overrep.results[,num.cols] = apply(overrep.results[,num.cols], 2, unlist)
		
		# Make plot of overrepresentation test p-values as the p-value for gene
		# inclusion changes. Only do for annotations with some p < 1
		
		# This should likely be spun out to its own function
		
		if (sum(overrep.results$overrep.p < 1, na.rm=TRUE) > 0) {
			
			if (dir == "lo") {
				p.cutoff.str = paste0(" (Low p-value, high ", 
										proper.names[[est.name]],
										")")
				p.cutoff.file.suffix = "lowP"
			} else {
				p.cutoff.str = paste0(" (High p-value, low ", 
										proper.names[[est.name]],
										")")
				p.cutoff.file.suffix = "highP"
			}
			
			pdf(paste0(out.dir, est.name, ".", this.anno, ".", 
						anno.val, ".", p.cutoff.file.suffix, ".pdf"))
				plot(overrep.results$overrep.p ~ overrep.results$p.val.cutoff,
					main=paste(proper.names[[this.anno]], anno.val, sep=" - "),
					xlab=paste0(proper.names[[est.name]], 
								" P-value cutoff for gene inclusion ",
								p.cutoff.str), 
					ylab="Overrepresentation test p-value",
					ylim=c(0,1),
					type='b', 
					pch=c(1, 16)[as.numeric(overrep.results$overrep.p < 0.05)+1])
				abline(h=0.05, lty=2)
			dev.off()
		}
		
		overrep.results.list[[dir]] = overrep.results
	}
	
	overrep.results.both = do.call(rbind, overrep.results.list)
	
	return(overrep.results.both)

}
