#!/usr/bin/env Rscript

# PANTHER parsing bit is based largely on code from this site:
# martinsbioblogg.wordpress.com/2013/02/11/using-r-accessing-panther-classifications/

# ----------------------------------------------------------------------------------------
# --- Download and parse PANTHER database for macaque
# ----------------------------------------------------------------------------------------

source("scripts/get.panther.db.macaque.R")

# ----------------------------------------------------------------------------------------
# --- Convert gene names to Ensembl IDs by parsing Uniprot reference proteome DB
# ----------------------------------------------------------------------------------------

source("scripts/parse.ref.proteome.R")

# ----------------------------------------------------------------------------------------
# --- Helper function to convert gene names to Ensembl IDs for enrichment input data
# ----------------------------------------------------------------------------------------

to.ensembl = function (est) {
	
	# uni will be the same as assoc.ids, but I am leaving it for sake of simplicity
	if (exists("uni") == FALSE) {
		uni = read.table("data/9544_macaca_mulatta.uniprot2genename")
		names(uni) = c("uniprot", "gene")
		uni$gene.cap = toupper(uni$gene)
	}
		
	est$gene.cap = toupper(est$gene)
		
	ests.idens = merge(x=uni, y=est, by.x="gene.cap", by.y="gene.cap")[,-3]
	
	names(ests.idens)[3] = "gene"

	return (ests.idens)
}

# ----------------------------------------------------------------------------------------
# --- Import function to query PANTHER data by Uniprot accession
# ----------------------------------------------------------------------------------------

source("scripts/query.panther.R")

# Test:
# print (query.panther(acc="F6QFW9", to.return="pathway"))

# ----------------------------------------------------------------------------------------
# --- Import function to annotate a gene list
# ----------------------------------------------------------------------------------------

source("scripts/annotate.gene.list.R")

# Test:
# print (annotate.gene.list(	gene.list=est.idens, 
#								uniprot.col.name="uniprot", 
#								which.anno="pathway"))

# ----------------------------------------------------------------------------------------
# --- Helper function to get matches given an annotation value
# ----------------------------------------------------------------------------------------

get.anno.matches = function (annotated.gene.list, anno.val) {

	matches = annotated.gene.list[which(annotated.gene.list$annotation == anno.val),]
	matches.nodup = matches[!duplicated(matches$gene.cap),]
	
	return(matches.nodup)
	
}

# ----------------------------------------------------------------------------------------
# --- Helper function to get list of all annotation values in PANTHER
# ----------------------------------------------------------------------------------------

get.all.panther.annos = function (this.anno, panther.classification) {
	
	if (this.anno == "pathway") {
		all.annotations = gsub("[#>]", "", 
					unique(sort(as.character(
							grep(">", unlist(panther.classification$panther.pathway), 
								value=TRUE)))))
	} else {

		all.annotations = gsub("^#", "", gsub("^GO:", "", 
					unique(sort(as.character(
							unlist(panther.classification[[this.anno]]))))))
	}
	
	return(all.annotations)
}

# ----------------------------------------------------------------------------------------
# --- Main function to classify genes and do enrichment or overrep tests
# ----------------------------------------------------------------------------------------

do.classify.and.test = function (	input.file,
									est.name, 
									output.prefix,
									this.anno=c("pathway", "go.mf", "go.bp", "go.cc", 
												"panther.proteinclass"),
									do.enrich.test=c(TRUE, FALSE),
									do.overrep.test=c(TRUE, FALSE)) {
	
	test.results = list()
	
	# ------------------------------------------------------------------------------------
	# --- Read in file for enrichment test
	# ------------------------------------------------------------------------------------
		
	if (do.enrich.test) {
	
		est = read.table(input.file, stringsAsFactors=FALSE)
	
		names(est) = c("gene", "estimate")
		
		est.idens = to.ensembl(est)
	}
	
	# ------------------------------------------------------------------------------------
	# --- Annotate gene list (enrichment test)
	# ------------------------------------------------------------------------------------

	if (do.enrich.test) {

		# Annotated gene list
		anno.est.val = annotate.gene.list(	gene.list=est.idens, 
											uniprot.col.name="uniprot", 
											which.anno=this.anno)
	}
	
	# Find all annotations in PANTHER
	all.annotations = get.all.panther.annos (this.anno, panther.classification) 
	
	# ------------------------------------------------------------------------------------
	# --- Do enrichment test
	# ------------------------------------------------------------------------------------

	# Bring in script that does enrichment test for one annotation value and all values
	source("scripts/enrich.test.annotation.R")
	source("scripts/enrich.test.all.R")
	
	# Do enrichment test for all values
	
	if (do.enrich.test) {
		anno.enrich.results = enrich.test.all (anno.est.val, all.annotations)
	}
	
	# ------------------------------------------------------------------------------------
	# --- Get names of annotations (pathway names, cell components, etc.)
	# ------------------------------------------------------------------------------------

	# Bring in necessary functions
	source("scripts/get.all.anno.vals.R")
	source("scripts/get.anno.name.R")
	
	anno.split = get.all.anno.vals(this.anno, panther)

	if (do.enrich.test) {

		anno.enrich.results$name = sapply(	rownames(anno.enrich.results), 
											get.anno.name,
											all.anno.vals = anno.split)
	}

	# ------------------------------------------------------------------------------------
	# --- Print out annotated gene list
	# ------------------------------------------------------------------------------------

	if (do.enrich.test) {
	
		anno.est.val.towrite = anno.est.val[!is.na(anno.est.val$annotation),]
		anno.est.val.towrite$name = sapply(	anno.est.val.towrite$annotation, 
											get.anno.name,
											all.anno.vals = anno.split)
		
		write.table(anno.est.val.towrite, paste0(output.prefix, ".annotated.list.txt"))
	}
	
	# ------------------------------------------------------------------------------------
	# --- Find significantly enriched annotation values
	# ------------------------------------------------------------------------------------
	
	# Get annotation values with p < 0.05

	if (do.enrich.test) {
		anno.sig = anno.enrich.results[which(anno.enrich.results$p < 0.05),]
	}
	
	# ------------------------------------------------------------------------------------
	# --- Plot cummulative density (ECDF) for significantly enriched values
	# ------------------------------------------------------------------------------------
	
	if (do.enrich.test) {

		# Should this be mini helper function? Also used in enrich.test.all.R
		all = anno.est.val[!duplicated(anno.est.val$gene.cap),]
	
		# This should be re-written as a separate function

		pdf(paste0(output.prefix, ".pdf"))
		
			# Plot overall ECDF
			all.ecdf  = ecdf(all$estimate)		
			color.counter = 2	
			plot(all.ecdf, 
					verticals = FALSE, do.points = TRUE,
					main=paste0('Empirical Cumluative Distribution - ', 
									proper.names[[est.name]], ' - ', 
									proper.names[[this.anno]]),
					xlab="Estimate", ylab="Fraction", cex=0.25)
		
			for (anno.val in row.names(anno.sig)) {
	
				matches.nodup = get.anno.matches(	annotated.gene.list = anno.est.val, 
													anno.val = anno.val)
				n.matches = nrow(matches.nodup)
				anno.ecdf = ecdf(matches.nodup$estimate)
				
				plot(anno.ecdf, add=TRUE, cex=0.5, col=color.counter)
				color.counter = color.counter + 1
			}
			
			if (nrow(anno.sig) > 0) {
				legend(x="bottomright", legend=anno.sig$name, col=2:(1+nrow(anno.sig)), 
						merge=FALSE, lty=1, pch=19, cex=0.75, bg="white")
			}
		
		dev.off()
	}
	
	# ------------------------------------------------------------------------------------
	# --- Print results to file
	# ------------------------------------------------------------------------------------
	
	if (do.enrich.test) {

		out.file = paste0(output.prefix, ".txt")
		write.table(anno.enrich.results, file=out.file)
	}
	
	# ------------------------------------------------------------------------------------
	# --- Store results in list
	# ------------------------------------------------------------------------------------
	
	if (do.enrich.test) {
		test.results[["enrich"]] = anno.enrich.results
	}
	
	# ------------------------------------------------------------------------------------
	# --- Read in file for overrep test
	# ------------------------------------------------------------------------------------
	
	# For the overrepresentation test, read in input file
	# containing gene names, number of windows gene occurs in, minimum estimate, 
	# maximum estimate, average estimate, and minimum p-value
	
	if (do.overrep.test) {
		extremes = read.table(input.file, stringsAsFactors=FALSE, header=TRUE)
	
		# Annotation just wants two columns, gene name and estimate
		extremes.est = extremes[,c(1,5)]
		names(extremes.est) = c("gene", "estimate")

		extremes.idens = to.ensembl(extremes.est)

	}
	
	# ------------------------------------------------------------------------------------
	# --- Annotate gene list (overrep test)
	# ------------------------------------------------------------------------------------

	if (do.overrep.test) {

		# Annotated gene list
		anno.extremes.val = annotate.gene.list(	gene.list=extremes.idens, 
												uniprot.col.name="uniprot", 
												which.anno=this.anno)
	}
	
	# ------------------------------------------------------------------------------------
	# --- Do overrep test
	# ------------------------------------------------------------------------------------

	# Bring in script that does overrep test for one annotation value,
	# and all p-value cutoffs for one annotation value
	source("scripts/overrep.test.annotation.R")
	source("scripts/overrep.test.annotation.pval.R")
		
	# Do overrep test for all values
	
	if (do.overrep.test) {
		anno.overrep.results = sapply(all.annotations, overrep.test.annotation,
										est.name=est.name, 
										annotated.gene.list=anno.extremes.val,
										extremes=extremes)
	}	

	# ------------------------------------------------------------------------------------
	# --- Write overrep test results to file
	# ------------------------------------------------------------------------------------

	# Bring in function to write overrep results to file
	source("scripts/write.overrep.test.results.R")
		
	if (do.overrep.test) {
		write.overrep.test.results(anno.overrep.results, est.name, this.anno, anno.split)
	}
	
	# ------------------------------------------------------------------------------------
	# --- Store results in list
	# ------------------------------------------------------------------------------------

	if (do.overrep.test) {
		test.results[["overrep"]] = anno.overrep.results
	}
	
	# ------------------------------------------------------------------------------------
	# --- Return enrichment and/or overrepresentation test results as list
	# ------------------------------------------------------------------------------------

	return(test.results)
}

# ----------------------------------------------------------------------------------------
# --- Call function on all input files (This part specific to Awash project.)
# ----------------------------------------------------------------------------------------

# Proper names for graphics
proper.names = list()

proper.names[["fst"]]   = "Fst"

proper.names[["fst.pval"]]   = "Fst P-value"

proper.names[["pathway"]] = "PANTHER Pathway"
proper.names[["go.mf"]]   = "GO Molecular Function"
proper.names[["go.bp"]]   = "GO Biological Process"
proper.names[["go.cc"]]   = "GO Cellular Component"
proper.names[["panther.proteinclass"]] = "PANTHER Protein Class"

# All possible annotation types
annotation.types = c("pathway", "go.mf", "go.bp", "go.cc", "panther.proteinclass")

# Output directory
out.dir = "results/panther_results/"
dir.create(file.path(out.dir), showWarnings = FALSE)

# Input files for enrichment test. Windowed:
input.enrich.win = Sys.glob("results/gene_lists/enrich_test/*.all.orthologs.genes.list.*")
input.enrich.win.est.names = gsub(".*enrich_test/(.*)\\.all.*", "\\1", input.enrich.win)

# Input files for enrichment test. Around ROI:
input.enrich.roi = Sys.glob("refGene/refGene.sort.gtf.papAnu2.*.sigma*.panther.*")
# Remove files that are results (output) of this script
input.enrich.roi = grep("panther.enrich.", input.enrich.roi, invert=TRUE, value=TRUE)
input.enrich.roi.est.names = gsub(".+papAnu2\\.([^\\.]+)\\.sigma.+", "\\1", 
									input.enrich.roi)

# Input files for overrep test. (Only) windowed:
input.overrep.win = Sys.glob(paste0("results/gene_lists/overrep_test/",
									"*.all.orthologs.genes.out.est.pval.extremes"))
input.overrep.win.est.names = gsub(".*overrep_test/(.+)\\.all.*", "\\1", 
									input.overrep.win)

# ----------------------------------------------------------------------------------------
# --- Call function to do enrichment test on results of windowed smoothing
# ----------------------------------------------------------------------------------------

for (est.idx in 1:length(input.enrich.win)) {

	input.est = input.enrich.win[est.idx]
	input.est.name = input.enrich.win.est.names[est.idx]
	
	for (anno.idx in 1:length(annotation.types)) {
	
		this.anno = annotation.types[anno.idx]
		
		results = do.classify.and.test (input.file = input.est,
										est.name = input.est.name, 
										output.prefix = paste0(out.dir, input.est.name, 
															".all.orthologs.genes.list.",
															"panther.enrich.",
															this.anno),
										this.anno = this.anno,
										do.enrich.test = TRUE,
										do.overrep.test = FALSE)
	}
}

# ----------------------------------------------------------------------------------------
# --- Call function to do enrichment test on results smoothed around regions of interest
# ----------------------------------------------------------------------------------------

for (est.idx in 1:length(input.enrich.roi)) {

	input.est = input.enrich.roi[est.idx]
	input.est.name = input.enrich.roi.est.names[est.idx]
	
	write(paste0("Processing ", input.est), stderr())
	
	for (anno.idx in 1:length(annotation.types)) {
	
		this.anno = annotation.types[anno.idx]

		# One unresolved problem: The overrep code ignores the output prefix.
		# For instance in the PDF generation in overrep.test.annotation.R
		# and the output text file generation in write.overrep.test.results.R
		results = do.classify.and.test (input.file = input.est,
										est.name = input.est.name, 
										output.prefix = paste0(input.est, 
																".panther.enrich.",
																this.anno),
										this.anno = this.anno,
										do.enrich.test = TRUE,
										do.overrep.test = FALSE)		
	}
}

# ----------------------------------------------------------------------------------------
# --- Call function to do overrep test on results of windowed smoothing
# ----------------------------------------------------------------------------------------

for (est.idx in 1:length(input.overrep.win)) {

	input.est = input.overrep.win[est.idx]
	input.est.name = input.overrep.win.est.names[est.idx]
	
	write(paste0("Processing ", input.est), stderr())
	
	for (anno.idx in 1:length(annotation.types)) {
	
		this.anno = annotation.types[anno.idx]

		results = do.classify.and.test (input.file = input.est,
										est.name = input.est.name, 
									output.prefix = paste0(input.est.name, 
															".overrep.panther.",
															this.anno),
									this.anno = this.anno,
									do.enrich.test = FALSE,
									do.overrep.test = TRUE)
	}
}
