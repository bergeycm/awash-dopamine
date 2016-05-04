#!/usr/bin/Rscript

library(xtable)
library(ggplot2)
library(reshape2)

options(scipen=999)

# Get directories for output ready. First windowed.
enrich.window.out.dir = "results/panther_results/parsed_results/windowed/enrich"
dir.create(enrich.window.out.dir, showWarnings=FALSE, recursive=TRUE)

overrep.out.dir = "results/panther_results/parsed_results/windowed/overrep"
dir.create(overrep.out.dir, showWarnings=FALSE, recursive=TRUE)

# The directories results/panther_results/parsed_results/windowed/enrich/ and
# results/panther_results/parsed_results/roi/enriched/ will contain:
# - panther_enrichment_outliers_*.tex
# And results/panther_results/parsed_results/windowed/overrep/ will contain
# - panther_overrep_outliers_*.tex

# Do same for KS around ROI enrichment results
enrich.roi.out.dir = "results/panther_results/parsed_results/roi/enrich"
dir.create(enrich.roi.out.dir, showWarnings=FALSE, recursive=TRUE)

# All tests
anno.types = c("go.bp", "go.cc", "go.mf", "pathway", "panther.proteinclass")

# Make human readable test names
anno.types.human = list()
anno.types.human[[anno.types[1]]] = "GO Biological Processes"
anno.types.human[[anno.types[2]]] = "GO Cellular Components"
anno.types.human[[anno.types[3]]] = "GO Molecular Functions"
anno.types.human[[anno.types[4]]] = "PANTHER Pathways"
anno.types.human[[anno.types[5]]] = "Protein Classes"

# All estimates
all.ests = c("fst")

# Make human readable estimate names
estimates.human = list()
estimates.human[[all.ests[1]]] = "Fst"

# ========================================================================================
# --- Read in PANTHER results files from statistical enrichment tests 
# ========================================================================================

# List to hold all PANTHER results
panther.out = list()

# List to hold PANTHER window enrichment test results
panther.out[["enrich_window"]] = list()

# List to hold PANTHER overrepresentation test results
panther.out[["overrep"]] = list()

# List to hold PANTHER ROI enrichment test results (the latter based on p-value)
panther.out[["enrich_ROI"]] = list()
panther.out[["enrich_ROI.pval"]] = list()

for (est.ind in 1:length(all.ests)) {

	this.est = all.ests[est.ind]
	
	panther.out[["enrich_window"]][[this.est]] = list()
	panther.out[["enrich_ROI"]][[this.est]] = list()
	panther.out[["enrich_ROI.pval"]][[this.est]] = list()
	
	# Create lists for each sigma for ROI enrichment test results
	for (sigma in c(50000, 100000, 150000)) {
		panther.out[["enrich_ROI"]][[this.est]][[sigma]]   = list()
		panther.out[["enrich_ROI.pval"]][[this.est]][[sigma]] = list()
	}
	
	for (test.ind in 1:length(anno.types)) {
	
		# --------------------------------------------------------------------------------
		# --- Enrichment test results - windowed
		# --------------------------------------------------------------------------------
				
		test = anno.types[test.ind]
		
		# PANTHER windowed enrichment test results in files like:
		# results/panther_results/fst.all.orthologs.genes.list.panther.enrich.go.bp.txt
		
		filename = paste0("results/panther_results/", this.est, 
								".all.orthologs.genes.list.panther.enrich.",
								test, ".txt")
		tmp = read.table(filename, header=TRUE)
		
		# Stupid Lipoate_biosynthesis has underscore and messes up the LaTeX compilation
		tmp$name = gsub("_", " ", tmp$name)
			
		# Put this data.frame in big results list
		panther.out[["enrich_window"]][[this.est]][[test]] = tmp
		
		rm("filename", "tmp")
		
		# --------------------------------------------------------------------------------
		# --- Enrichment test results - ROI
		# --------------------------------------------------------------------------------
						
		# PANTHER ROI enrichment test results in files like:
		# refGene/refGene.sort.gtf.papAnu2.fst.sigma100000.txt
		
		for (sigma in c(50000, 100000, 150000)) {
				
			filename = paste0("refGene/refGene.sort.gtf.papAnu2.", this.est,
									".sigma", sigma, ".panther.", this.est, 
									".panther.enrich.", test, ".txt")
			tmp = read.table(filename, header=TRUE)
			
			tmp$name = gsub("_", " ", tmp$name)

			# Put this data.frame in big results list
			panther.out[["enrich_ROI"]][[this.est]][[test]][[sigma]] = tmp
			
			# Also parse PANTHER results based on p-value
			filename = paste0("refGene/refGene.sort.gtf.papAnu2.", this.est,
									".sigma", sigma, ".panther.", this.est, ".pval",
									".panther.enrich.", test, ".txt")
			
			tmp = read.table(filename, header=TRUE)
			
			tmp$name = gsub("_", " ", tmp$name)
			
			# Put this data.frame in big results list
			panther.out[["enrich_ROI.pval"]][[this.est]][[test]][[sigma]] = tmp
			
			rm("filename", "tmp")
		}
		
		# --------------------------------------------------------------------------------
		# --- Overrepresentation test results
		# --------------------------------------------------------------------------------
		
		# PANTHER overrep test results in files like:
		# results/panther_results/fst.overrep.go.bp.txt

		filename = paste0("results/panther_results/", this.est, 
								".overrep.",
								test, ".txt")
		tmp = read.table(filename, header=TRUE, stringsAsFactor=FALSE)
		
		tmp$name = gsub("_", " ", tmp$name)

		# Put this data.frame in big results list
		panther.out[["overrep"]][[this.est]][[test]] = tmp
	
		rm("filename", "tmp")
	}
}

# ========================================================================================
# --- Output tables of enrichment and overrep test outliers
# ========================================================================================

for (est.ind in 1:length(all.ests)) {

	this.est = all.ests[est.ind]
	
	for (test.ind in 1:length(anno.types)) {
	
		test = anno.types[test.ind]
		
		# --------------------------------------------------------------------------------
		# --- Output enrichment test results in LaTeX table - windowed
		# --------------------------------------------------------------------------------		
	
		curr.results = panther.out[["enrich_window"]][[this.est]][[test]]
		
		# Get significant results (p < 0.05)
		col.order = c(5,1,3,4)
		enrich.sig.5 = curr.results[which(curr.results$p < 0.05), col.order]
		
		# Sort by p-value
		enrich.sig.5 = enrich.sig.5[order(enrich.sig.5$p),]
		
		# Round p-value to three significant digits
		enrich.sig.5$p = signif(enrich.sig.5$p, 3)
		
		# - Make column headers human readable -------------------------------------------
	
		names(enrich.sig.5) = 	c(anno.types.human[[test]], 
								"Gene Count", "+/$-$", "Enrichment $P$")
	
		# - Change it so that "-" is in Math format: "$-$" -------------------------------
		
		levels(enrich.sig.5[,3]) = paste("$", levels(enrich.sig.5[,3]), "$", sep="")
		
		# - Put two stars on values < 0.01 and a star on values < 0.05 -------------------
		
		one.star = (enrich.sig.5[,4] <= 0.05) & (enrich.sig.5[,4] > 0.01)
		two.star = enrich.sig.5[,4] <= 0.01	
		
		enrich.sig.5[,4][one.star] = paste("*",  enrich.sig.5[,4][one.star], sep="")
		enrich.sig.5[,4][two.star] = paste("**", enrich.sig.5[,4][two.star], sep="")
	
		# Write LaTeX file ---------------------------------------------------------------
	
		tex.out = paste0(enrich.window.out.dir, "/panther_enrichment_outliers_window_", 
							this.est, "_", test, ".tex")
		xt = xtable(enrich.sig.5)
		align(xt) = c('l', 'p{3in}||', rep('r', 3))
			
		sink(tex.out)
		
			cat("\\documentclass{article}", 
				"\\usepackage{graphicx}",
				"\\DeclareGraphicsExtensions{.pdf}",
				"\\usepackage[top=0.5in, bottom=0.5in, left=0.5in, right=0.5in]{geometry}",
				"\\begin{document}", sep="\n")
			
			print.xtable(xt, 
						include.rownames = FALSE, 
						sanitize.text.function = function(x){x}, 
						size="scriptsize")
			
			# Add graph from PANTHER
			graph.path = paste0("../../../", this.est, 
								".all.orthologs.genes.list.panther.enrich.", test)
			
			cat(paste0("\\centerline{",
						"\\includegraphics[width=0.75\\textwidth]",
						"{{", graph.path, "}.pdf}",
						"}"), sep="\n")
	
	
			cat("\\end{document}", sep="\n")
			
		sink()
		
		# --------------------------------------------------------------------------------
		# --- Output enrichment test results in LaTeX table - ROI
		# --------------------------------------------------------------------------------		
	
		# Do twice: once for weighted estimates, one for p-values
		for (p.val.var in c("", ".pval")) {
		
			for (sigma in c(50000, 100000, 150000)) {
				
				curr.results = panther.out[[paste0("enrich_ROI", p.val.var)]][[this.est]][[test]][[sigma]]
				
				# Get significant results (p < 0.05)
				col.order = c(5,1,3,4)
				enrich.sig.5 = curr.results[which(curr.results$p < 0.05), col.order]
			
				# Sort by p-value
				enrich.sig.5 = enrich.sig.5[order(enrich.sig.5$p),]
			
				# Round p-value to three significant digits
				enrich.sig.5$p = signif(enrich.sig.5$p, 3)
			
				# - Make column headers human readable -------------------------------------------
				
				names(enrich.sig.5) = 	c(anno.types.human[[test]], 
									"Gene Count", "+/$-$", "Enrichment $P$")
		
				# - Change it so that "-" is in Math format: "$-$" -------------------------------
			
				levels(enrich.sig.5[,3]) = paste("$", levels(enrich.sig.5[,3]), "$", sep="")
			
				# - Put two stars on values < 0.01 and a star on values < 0.05 -------------------
			
				one.star = (enrich.sig.5[,4] <= 0.05) & (enrich.sig.5[,4] > 0.01)
				two.star = enrich.sig.5[,4] <= 0.01	
			
				enrich.sig.5[,4][one.star] = paste("*",  enrich.sig.5[,4][one.star], sep="")
				enrich.sig.5[,4][two.star] = paste("**", enrich.sig.5[,4][two.star], sep="")
		
				# Write LaTeX file ---------------------------------------------------------------
		
				tex.out = paste0(enrich.roi.out.dir, "/panther_enrichment_outliers_ROI_", 
									this.est, p.val.var, "_", test, "_", sigma, ".tex")
				xt = xtable(enrich.sig.5)
				align(xt) = c('l', 'p{3in}||', rep('r', 3))
				
				sink(tex.out)
			
					cat("\\documentclass{article}", 
						"\\usepackage{graphicx}",
						"\\DeclareGraphicsExtensions{.pdf}",
						"\\usepackage[top=0.5in, bottom=0.5in, left=0.5in, right=0.5in]{geometry}",
						"\\begin{document}", sep="\n")
					
					print.xtable(xt, 
								include.rownames = FALSE, 
								sanitize.text.function = function(x){x}, 
								size="scriptsize")
					
					# Add graph from PANTHER
					graph.path = paste0("../../../../../refGene/",
										"refGene.sort.gtf.papAnu2.", this.est,
										".sigma", sigma, ".panther.", this.est, p.val.var,
										".panther.enrich.", test)
					
					cat(paste0("\\centerline{",
								"\\includegraphics[width=0.75\\textwidth]",
								"{{", graph.path, "}.pdf}",
								"}"), sep="\n")
			
			
					cat("\\end{document}", sep="\n")
					
				sink()
			}
		}
		
		# --------------------------------------------------------------------------------
		# --- Output enrichment test results in LaTeX table - ROI (all sigmas together)
		# --------------------------------------------------------------------------------		
	
		# Do twice: once for weighted estimates, one for p-values
		for (p.val.var in c("", ".pval")) {
		
			get.this.sigma.results = function(sigma) {
				return(cbind(panther.out[[paste0("enrich_ROI", p.val.var)]]
										[[this.est]][[test]][[sigma]], sigma))
			}
		
			curr.results = do.call(rbind, 
							lapply(c(50000, 100000, 150000), 
								function(x) {get.this.sigma.results(x)}))
				
			# Get significant results (p < 0.05)
			col.order = c(5,6,1,3,4)
			enrich.sig.5 = curr.results[which(curr.results$p < 0.05), col.order]
			
			# Figure out minimum p-value for each annotation for sorting
			if (nrow(enrich.sig.5) == 0) {
				next
			}
			p.mins = aggregate(p ~ name, enrich.sig.5, min)
			names(p.mins)[2] = "min.p"
			enrich.sig.5 = merge(enrich.sig.5, p.mins, by="name")

			# Sort by minimum p-value, then name, then sigma
			enrich.sig.5 = enrich.sig.5[order(	enrich.sig.5$min.p, 
												enrich.sig.5$name, 
												enrich.sig.5$sigma),]
			# Remove min.p column
			enrich.sig.5 = enrich.sig.5[,-c(6)]
			
			# Round p-value to three significant digits
			enrich.sig.5$p = signif(enrich.sig.5$p, 3)
			
			# - Make column headers human readable ---------------------------------------
				
			names(enrich.sig.5) = 	c(anno.types.human[[test]], "$\\sigma$",
									"Gene Count", "+/$-$", "Enrichment $P$")
		
			# - Change it so that "-" is in Math format: "$-$" ---------------------------
			
			levels(enrich.sig.5[,4]) = paste("$", levels(enrich.sig.5[,4]), "$", sep="")
			
			# - Put two stars on values < 0.01 and a star on values < 0.05 ---------------
			
			one.star = (enrich.sig.5[,5] <= 0.05) & (enrich.sig.5[,5] > 0.01)
			two.star = enrich.sig.5[,5] <= 0.01	
			
			enrich.sig.5[,5][one.star] = paste("*",  enrich.sig.5[,5][one.star], sep="")
			enrich.sig.5[,5][two.star] = paste("**", enrich.sig.5[,5][two.star], sep="")
			
			# - Deal with duplicated annotations -----------------------------------------

			# Figure out last lines for each annotation, to draw lines between then
			breaks = which(enrich.sig.5[,1] != c(enrich.sig.5[c(-1),1], NA))
		
			# Figure out which annotations do not match the line above, 
			# to keep those and delete redundant labels
			firsts = c(1,which(c(NA, enrich.sig.5[,1]) != c(enrich.sig.5[,1], NA)))
		
			# Delete redundant labels and counts of genes in reference list
			# Replace with ditto mark		
			if (nrow(enrich.sig.5) > 0) {
				enrich.sig.5[-firsts,1] = '"'
			}
		
			# - Write LaTeX file ---------------------------------------------------------
		
			tex.out = paste0(enrich.roi.out.dir, "/panther_enrichment_outliers_ROI_", 
								this.est, p.val.var, "_", test, "_allsigmas.tex")
			xt = xtable(enrich.sig.5, digits = c(0,0,0,0,0,NA))
			align(xt) = c('l', 'p{3in}||', rep('r', 4))
			
			sink(tex.out)
			
				cat("\\documentclass{article}", 
					"\\usepackage{graphicx}",
					"\\DeclareGraphicsExtensions{.pdf}",
					"\\usepackage[top=0.5in, bottom=0.5in, left=0.5in, right=0.5in]{geometry}",
					"\\begin{document}", sep="\n")
				
				print.xtable(xt, 
							include.rownames = FALSE, 
							sanitize.text.function = function(x){x}, 
							size="scriptsize")
				
				cat("\\end{document}", sep="\n")
				
			sink()
		}
		
		# --------------------------------------------------------------------------------
		# --- Output overrep test results in LaTeX table
		# --------------------------------------------------------------------------------

		curr.results = panther.out[["overrep"]][[this.est]][[test]]
		
		# Get significant results (p < 0.05)
		### CHECK
		col.order = c(12,11,4:10)
		overrep.sig.5 = curr.results[which(curr.results$overrep.p < 0.05), col.order]

		# Sort by annotation name
		### CHECK
		overrep.sig.5 = overrep.sig.5[order(overrep.sig.5$name),]
		
		# Round p-value to three significant digits
		overrep.sig.5$overrep.p = signif(overrep.sig.5$overrep.p, 3)
		
		# Split into high and low values
		overrep.sig.5.hi = overrep.sig.5[which(overrep.sig.5$hi.or.lo.p.cutoff == 'hi'),]
		overrep.sig.5.lo = overrep.sig.5[which(overrep.sig.5$hi.or.lo.p.cutoff == 'lo'),]		
		
		# Remove column with hi or lo
		overrep.sig.5.hi = overrep.sig.5.hi[,-c(2)]
		overrep.sig.5.lo = overrep.sig.5.lo[,-c(2)]
		
		# - Make column headers human readable -------------------------------------------
	
		overrep.names = c(anno.types.human[[test]], 
							paste0(estimates.human[[this.est]], " ($\\alpha$ cutoff)"), 
							"Obs (Ref.)", "Obs", "Exp", 
							"+/$-$", "Fold Enrich", "Overrep $P$")
	
		names(overrep.sig.5.hi) = names(overrep.sig.5.lo) = overrep.names
		
		# High p-values (Low estimates) should have "Low" before estimate name
		# and vice versa for Low p-values (High estimates)
		names(overrep.sig.5.hi)[2] = paste("Low",  names(overrep.sig.5.hi)[2])
		names(overrep.sig.5.lo)[2] = paste("High", names(overrep.sig.5.lo)[2])
	
		# - Change it so that "-" is in Math format: "$-$" -------------------------------
		
		if (nrow(overrep.sig.5.hi) > 0) {
			overrep.sig.5.hi[,6] = paste0("$", overrep.sig.5.hi[,6], "$")
		}
		if (nrow(overrep.sig.5.lo) > 0) {
			overrep.sig.5.lo[,6] = paste0("$", overrep.sig.5.lo[,6], "$")
		}
		
		# - Put two stars on values < 0.01 and a star on values < 0.05 -------------------
		
		one.star.hi = (overrep.sig.5.hi[,8] <= 0.05) & (overrep.sig.5.hi[,8] > 0.01)
		one.star.lo = (overrep.sig.5.lo[,8] <= 0.05) & (overrep.sig.5.lo[,8] > 0.01)
		two.star.hi = overrep.sig.5.hi[,8] <= 0.01	
		two.star.lo = overrep.sig.5.lo[,8] <= 0.01	
		
		overrep.sig.5.hi[,8][one.star.hi] = paste0("*",  overrep.sig.5.hi[,8][one.star.hi])
		overrep.sig.5.lo[,8][one.star.lo] = paste0("*",  overrep.sig.5.lo[,8][one.star.lo])
		
		overrep.sig.5.hi[,8][two.star.hi] = paste0("**", overrep.sig.5.hi[,8][two.star.hi])
		overrep.sig.5.lo[,8][two.star.lo] = paste0("**", overrep.sig.5.lo[,8][two.star.lo])
		
		# Figure out last lines for each annotation, to draw lines between then
		breaks.hi = which(overrep.sig.5.hi[,1] != c(overrep.sig.5.hi[c(-1),1], NA))
		breaks.lo = which(overrep.sig.5.lo[,1] != c(overrep.sig.5.lo[c(-1),1], NA))
		
		# Figure out which annotations do not match the line above, 
		# to keep those and delete redundant labels
		firsts.hi = c(1,which(c(NA, overrep.sig.5.hi[,1]) != c(overrep.sig.5.hi[,1], NA)))
		firsts.lo = c(1,which(c(NA, overrep.sig.5.lo[,1]) != c(overrep.sig.5.lo[,1], NA)))
		
		# Delete redundant labels and counts of genes in reference list
		# Replace with ditto mark		
		if (nrow(overrep.sig.5.hi) > 0) {
			overrep.sig.5.hi[-firsts.hi,c(1,3)] = '"'
		}
		if (nrow(overrep.sig.5.lo) > 0) {
			overrep.sig.5.lo[-firsts.lo,c(1,3)] = '"'
		}
				
		# Write LaTeX file ---------------------------------------------------------------
	
		display = c('s','s','f','s','s','f','s','f','s')
		digits  = c( 0,  0,  3,  0,  0,  2,  0,  2,  0)
		
		tex.out = paste0(overrep.out.dir, "/panther_overrep_outliers_", 
							this.est, "_", test, ".tex")
		xt.hi = xtable(overrep.sig.5.hi, 
						caption=paste0("Significantly overrepresented ", 
										anno.types.human[[test]], 
										" categories found for high estimates of ", 
										estimates.human[[this.est]],
										"."),
						display=display, digits=digits)
		
		xt.lo = xtable(overrep.sig.5.lo, 
						caption=paste0("Significantly overrepresented ", 
										anno.types.human[[test]], 
										" categories found for low estimates of ", 
										estimates.human[[this.est]],
										"."),
						display=display, digits=digits)
						
		align(xt.hi) = align(xt.lo) = c('l', 'p{2in}||', rep('r', 7))
		
		addtorow     = list()
		addtorow$pos = list()
		addtorow$pos[[1]] = c(0)
		addtorow$command = c(paste0("\\hline \n",
									"\\endhead \n",
									"\\hline \n",
									"{\\footnotesize ",
										"(Table continued on next page)",
									"} \n",
									"\\endfoot \n",
									"\\endlastfoot \n"))		
		sink(tex.out)
		
			cat("\\documentclass{article}", 
				"\\usepackage{graphicx}",
				"\\usepackage{longtable}",
				"\\DeclareGraphicsExtensions{.pdf}",
				"\\usepackage[top=0.5in, bottom=0.5in, left=0.5in, right=0.5in]{geometry}",
				"\\usepackage{caption}",
				"\\captionsetup[table]{labelformat=empty}",
				"\\begin{document}", sep="\n")
			
			if (nrow(xt.hi) > 0) {
			
				print.xtable(xt.hi, 
							include.rownames = FALSE, 
							sanitize.text.function = function(x){x}, 
							size="scriptsize",
							tabular.environment = 'longtable', floating = FALSE,
							add.to.row = addtorow, 
							hline.after = c(-1, breaks.hi),
							caption.placement = "top")
			} else {
				cat (paste0("\\par\\scriptsize{",
							"(No significantly overrepresented ", 
							anno.types.human[[test]], 
							" categories found for high estimates of ", 
							estimates.human[[this.est]],
							".)}"), sep="\n")
			}		
			cat("\\", sep="\n")

			if (nrow(xt.lo) > 0) {

				print.xtable(xt.lo, 
							include.rownames = FALSE, 
							sanitize.text.function = function(x){x}, 
							size="scriptsize",
							tabular.environment='longtable', floating=FALSE,
							add.to.row = addtorow, 
							hline.after = c(-1, breaks.lo),
							caption.placement = "top")
			} else {
				cat (paste0("\\par\\scriptsize{",
							"(No significantly overrepresented ", 
							anno.types.human[[test]], 
							" categories found for low estimates of ", 
							estimates.human[[this.est]],
							".)}"), sep="\n")
			}
			cat("\\end{document}", sep="\n")
			
		sink()
	
	}
}

# Compile *.tex files outside of R

# cd results/panther_results/parsed_results/windowed/enrich/
# for tex in $( ls panther_enrichment_outliers_*.tex ); do
#     pdflatex $tex
# done
# cd ../../../../..

# cd results/panther_results/parsed_results/windowed/overrep/
# for tex in $( ls panther_overrep_outliers_*.tex ); do
#     pdflatex $tex
# done
# cd ../../../../..

# cd results/panther_results/parsed_results/roi/enrich/
# for tex in $( ls panther_enrichment_outliers_*.tex ); do
#     pdflatex $tex
# done
# cd ../../../../..

