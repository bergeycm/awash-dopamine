#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------------
# --- Do overrepresentation test for one annotation value and one p-value cutoff
# ------------------------------------------------------------------------------------

overrep.test.annotation.pval = function (p.lim, high.or.low.p=c("hi", "lo"), 
											extremes, all, matches.nodup, 
											pC, nC) {
		
	# Get test list (outliers) for this p-value cutoff
	# Depending on mode, this will grab either the low or high p-values
	# - Low  p-value indicates high estimate (Fst, etc.)
	# - High p-value indicates low estimate  (Fst, etc.)
	if (high.or.low.p == "lo") {
		all.outliers = extremes[extremes$Min.p < p.lim,]$Gene
	} else if (high.or.low.p == "hi") {
		all.outliers = extremes[extremes$Min.p > (1-p.lim),]$Gene
	}
	
	# Reduce to just those that are in PANTHER
	outliers = all.outliers[toupper(all.outliers) %in% all$gene.cap]
	
	# Test list (outliers) contains how many genes
	K = length(outliers)

	# Expected number of genes involved in category
	# Expectation:
	# The probability of observing a gene w/ a particular annotation
	# in the list of outliers is the same as in the reference list. 
	expect = K * pC

	# Actual observed genes in category
	kC = sum(toupper(outliers) %in% matches.nodup$gene.cap)
	
	# Over or underrepresented
	over.under = "";
	if (kC > expect) {
		over.under = "+"
	} else if (kC < expect) {
		over.under = "-"
	}
	
	# Fold enrichment
	fold.enrich = kC / expect

	# Do binomial test to get p-value
	# Make sure that K is positive and greater than kC
	if (K <= 0 || K <= kC) {
		overrep.p = NA
	} else {
		overrep.p = binom.test(kC, K, pC, alternative = "g")$p.value
	}
	
	return(list(p.lim, nC, kC, expect, over.under, fold.enrich, overrep.p))
}
