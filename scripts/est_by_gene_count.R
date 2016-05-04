#!/usr/bin/Rscript

options(stringsAsFactors = FALSE)
library(ggplot2)
library(plyr)

analyses = sub(	"results/(.*)\\.all.+", 
				"\\1", 
				Sys.glob("results/*.all.orthologs.genes.out"), 
				perl=TRUE)

for(analysis in analyses) {

	est.genes = read.table(	paste0("results/", analysis, ".all.orthologs.genes.out"), 
							sep="\t")
	names(est.genes) = c("baboon", "rhesus", "est", "pval", "genes")
	
	genes = strsplit(est.genes$genes, ";", fixed=TRUE)
	
	est.genes$gene.count = do.call(rbind, lapply(genes, function(x) length(x)))
	
	# Figure out which expression to use in plots
	if (analysis == "fst") {
		est.expr = expression(F[ST])
	} else {
		est.expr = "Kernel Smoothed Estimate"
	}
	
	# ========================================================================================
	# --- Linear regression
	# ========================================================================================
	
	gene.lm = lm(est ~ gene.count, data=est.genes)
	step.change = gene.lm$coefficients[["gene.count"]]
	
	sink(file=paste0("reports/", analysis, "_by_gene_count_lm.txt"))
		cat(paste("--- Linear regression with all gene counts: ---"), sep="\n")
		summary(gene.lm)
		cat(paste("Estimate change with each additional gene:", step.change), sep="\n")
	sink()
	
	# Get max for later use in ggplot (hard to get later after factoring)
	gene.ct.max = max(est.genes$gene.count)
	
	# --- Re-do with just regions that have gene counts from 0 to 4
	
	gene.lm.lt5 = lm(est ~ gene.count, data=est.genes[est.genes$gene.count < 5,])
	step.change.lt5 = gene.lm.lt5$coefficients[["gene.count"]]
	
	sink(file=paste0("reports/", analysis, "_by_gene_count_lm.txt"), append=TRUE)
		cat(paste("--- Linear regression with gene counts less than 5: ---"), sep="\n")
		summary(gene.lm.lt5)
		cat(paste("Estimate change with each additional gene:", step.change.lt5), sep="\n")
	sink()
	
	# ========================================================================================
	# --- Get average estimate by gene count
	# ========================================================================================
	
	est.summary = ddply(est.genes[is.na(est.genes$est) == FALSE,], ~gene.count,
						summarize, mean=mean(est), sd=sd(est))
	write.table(est.summary, file=paste0("reports/", analysis, "_by_gene_count_avg.txt"))
	
	# ========================================================================================
	# --- Create box plots with lm lines
	# ========================================================================================
	
	# --- Factoring for ggplot
	
	# Factor to keep levels without data from being excluded
	est.genes$gene.count = factor(est.genes$gene.count, 
									levels=seq(from=0, to=max(est.genes$gene.count)))
	
	# --- Plot all data
	
	p = ggplot(est.genes, aes(gene.count, est)) +
		geom_boxplot() + 
		stat_smooth(method="lm", se=FALSE, aes(group=1)) +
		scale_x_discrete("Genes in region", breaks=factor(0:gene.ct.max), drop=FALSE) +
		labs(y=est.expr)
	
	ggsave(filename=paste0("results/", analysis, "_by_gene_count_boxplot_lm.pdf"), 
			plot=p)
	
	# --- Re-do with just regions that have gene counts from 0 to 4
	
	q = ggplot(est.genes[est.genes$gene.count %in% 0:4,], aes(gene.count, est)) +
		geom_boxplot() + 
		stat_smooth(method="lm", se=FALSE, aes(group=1)) +
		scale_x_discrete("Genes in region (Truncated to regions with 0-4 genes)", 
			breaks=factor(0:gene.ct.max), drop=TRUE) +
		labs(y=est.expr)
	
	ggsave(filename=paste0("results/", analysis, "_by_gene_count_boxplot_lm_lt5.pdf"), 
			plot=q)
}
