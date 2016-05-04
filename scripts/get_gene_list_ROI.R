#!/usr/bin/env Rscript

# ========================================================================================
# --- Script to take results of kernel smoothing around ROI (such as refGenes) and
# --- create input file for PANTHER.
# ========================================================================================

# Usage:  Rscript scripts/get_gene_list_ROI.R [ROI.bed] [out.label]
# Ex:     Rscript scripts/get_gene_list_ROI.R "refGene/refGene.sort.gtf.papAnu2.bed" "fst"
# Input:  ROI BED file 
#			(e.g. refGene/refGene.sort.gtf.papAnu2.bed)
#         Kernel smoothing around ROI results files
#			(e.g. refGene/refGene.sort.gtf.papAnu2.fst.sigma150000.txt)
# Output: File for PANTHER enrichment test containing Fst values
#			(e.g. refGene/refGene.sort.gtf.papAnu2.fst.sigma150000.panther.fst)
#         File for PANTHER enrichment test containing p values
#			(e.g. refGene/refGene.sort.gtf.papAnu2.fst.sigma150000.panther.pval)

# ----------------------------------------------------------------------------------------
# --- Get user input
# ----------------------------------------------------------------------------------------

args = commandArgs(trailingOnly = TRUE)

# BED file of regions of interest
roi.bed.file = args[1]			# "refGene/refGene.sort.gtf.papAnu2.bed"

# Label used in naming of smoothed output files (e.g. "fst")
out.label = args[2]				# "fst"

# ----------------------------------------------------------------------------------------
# --- Read in input files
# ----------------------------------------------------------------------------------------

# Read in original BED file
roi.bed = read.table(roi.bed.file, header=TRUE)

ks.ests = list.files(path="refGene/", pattern=paste0(out.label, ".sigma\\d+.txt$"))

sigmas = gsub(".*sigma(\\d+).*", "\\1", ks.ests)

# Read in estimates that have been kernel smoothed around the ROIs

for (sigma in sigmas) {
	
	ks.in = paste0(	gsub(".bed", paste0(".", out.label), roi.bed.file), 
					".sigma", sigma, ".txt")
	this.est.ks.roi = read.table(ks.in, header=TRUE)
	names(this.est.ks.roi)[4:5] = paste0("sigma", sigma, ".", names(this.est.ks.roi)[4:5])
	
	roi.bed = merge(x=roi.bed, y=this.est.ks.roi, 
					by.x=c("Chr", "Begin", "End"), by.y=c("chr", "start", "end"))

}

# ----------------------------------------------------------------------------------------
# --- Create plots to see concordance between runs with different sigma values
# ----------------------------------------------------------------------------------------

ests.cols = grep(".ests",  names(roi.bed))
pval.cols = grep(".p.val", names(roi.bed))

# Plot smoothed values of runs with different sigma values
plot(roi.bed[,ests.cols])

# Plot p-values of runs with different sigma values
plot(roi.bed[,pval.cols])

# ----------------------------------------------------------------------------------------
# --- Output files for PANTHER enrichment analysis
# ----------------------------------------------------------------------------------------
 
# Output two files per sigma value, one with smoothed estimate and one with p-value

gene.col = which(names(roi.bed) == "GeneID")

for (sigma in sigmas) {
	
	panther = paste0(	gsub(".bed", paste0(".", out.label), roi.bed.file),
						".sigma", sigma, ".panther")
	panther.ests = paste0(panther, ".", out.label)
	panther.pval = paste0(panther, ".", out.label, ".pval")

	ests.col = grep(paste0("sigma", sigma, ".ests.wavg"), names(roi.bed))
	pval.col = grep(paste0("sigma", sigma, ".p.val"),     names(roi.bed))
	
	write.table(roi.bed[is.na(roi.bed[,ests.col]) == FALSE, c(gene.col, ests.col)], 
				file=panther.ests, 
				quote=FALSE, sep="\t", 
				row.names=FALSE, col.names=FALSE)

	write.table(roi.bed[is.na(roi.bed[,pval.col]) == FALSE, c(gene.col, pval.col)], 
				file=panther.pval, 
				quote=FALSE, sep="\t", 
				row.names=FALSE, col.names=FALSE)
}
