#!/usr/bin/Rscript

options(stringsAsFactors = FALSE)

library(grid)
library(scales)
library(xtable)

# Which expression to use in plots
est.expr = list()
est.expr[["fst"]]   = expression("Smoothed"~F[ST])

# ========================================================================================
# --- Plot window smoothed values across genome
# ========================================================================================

# ----------------------------------------------------------------------------------------
# --- Get chromosome lengths
# ----------------------------------------------------------------------------------------

chrom.lens = read.table("genomes/papAnu2/papAnu2.fai")
chrom.lens = chrom.lens[,1:2]
names(chrom.lens) = c("chr", "len")
# Max width of plots, in inches
max.size = 6.5
chrom.lens$scaled = max.size * (chrom.lens$len / max(chrom.lens$len))

# ----------------------------------------------------------------------------------------
# --- Read in gene locations
# ----------------------------------------------------------------------------------------

genes = list()
for (chr in 1:20) {
	gene.in.file = paste0("results/orthologs/orthologs.chr", chr, ".genes.out")
	this.gene = read.table(gene.in.file, header=FALSE, sep="\t")
	names(this.gene) = c("region.pap", "region.mac", "gene.list")
	
	# Get middle of window
	this.gene$middle = rowMeans(apply(do.call(rbind, 
				strsplit(this.gene$region.pap, split="[:-]"))[,2:3], 2, as.numeric))
	# Get gene count
	this.gene$gene.count = lapply(strsplit(this.gene$gene.list, ";"), length)

	genes[[chr]] = this.gene
}

# ----------------------------------------------------------------------------------------
# --- Read in window smoothed data for Fst
# ----------------------------------------------------------------------------------------

ests = c('fst')

win.smooth = list()
for (est in ests) {
	win.smooth[[est]]   = list()
	
	for(chr in 1:20) {
		est.in.file = paste0("results/", est, "_weighted_averages_chr", chr, ".txt")
		win.smooth[[est]][[chr]] = read.table(est.in.file, header=TRUE)
	}
}

# ----------------------------------------------------------------------------------------
# --- Get range of values for each estimate
# ----------------------------------------------------------------------------------------

ylims = list()
ylims[["fst"]]   = c(0,1)

# ----------------------------------------------------------------------------------------
# --- Plot window smoothed data for Fst across the genome
# ----------------------------------------------------------------------------------------

for (est in ests) {

	pdf(file=paste0("results/smoothed_", est, "_against_genome.pdf"),
		height=9, width=6.5)
	
	# Changing width fails here.
	layout(matrix(1:6), widths=chrom.lens$scaled, heights=rep(2, 6))
	par(mar=c(4, 4, 1, 1) + 0.1)
	
	for (chr in 1:20) {

		curr = win.smooth[[est]][[chr]]
		
		plot(curr$w.avg ~ curr$pos, type='n', 
			xlab=paste("Chromosome", chr), ylab=est.expr[[est]],
			ylim=ylims[[est]], xlim=c(0, max(win.smooth[[est]][[1]]$pos)))
		
		abline(v=curr$pos[curr$p.val < 0.05], col=alpha("yellow",0.3))
		abline(v=curr$pos[curr$p.val > 0.95], col=alpha("orange",0.3))
		points(curr$w.avg ~ curr$pos, type='l')
		
		# Draw over p-value highlight lines
		box(); box(); box()
		
		# Add genes
		gene.height = ylims[[est]][1] + (0.8 * (ylims[[est]][2] - ylims[[est]][1]))
		points(x = genes[[chr]]$middle, y=rep(gene.height, nrow(genes[[chr]])), 
			cex=unlist(genes[[chr]]$gene.count) / 4,
			pch=16)
		
	}
	
	dev.off()
}	
		
# ========================================================================================
# --- Plot both window smoothed and ROI smoothed values across genome
# ========================================================================================

# ----------------------------------------------------------------------------------------
# --- Read in window smoothed data for Fst
# ----------------------------------------------------------------------------------------

sigmas = c('50000', '100000', '150000')

roi.smooth = list()

for (est in ests) {

	roi.smooth[[est]] = list()
		
	for (sigma in sigmas) {
		
		est.in.file = paste0("refGene/refGene.sort.gtf.papAnu2.", est, 
									".sigma", sigma, ".txt")
		roi.smooth[[est]][[sigma]] = read.table(est.in.file, header=TRUE)
		roi.smooth[[est]][[sigma]]$sigma = sigma
		
		# Get middle of gene
		roi.smooth[[est]][[sigma]]$middle = (roi.smooth[[est]][[sigma]]$start + 
			roi.smooth[[est]][[sigma]]$end) / 2
		
		# Sort genes by location
		roi.smooth[[est]][[sigma]] = roi.smooth[[est]][[sigma]][
			order(roi.smooth[[est]][[sigma]]$chr, roi.smooth[[est]][[sigma]]$middle),]
	}
	
	roi.smooth[[est]] = do.call(rbind, roi.smooth[[est]])
}

# ----------------------------------------------------------------------------------------
# --- Plot both window smoothed data and ROI smoothed data (multiple sigmas)
# ----------------------------------------------------------------------------------------

for (est in ests) {

	pdf(file=paste0("results/smoothed_", est, "_against_genome_window_vs_roi.pdf"),
		height=9, width=6.5)
	
	# Changing width fails here.
	layout(matrix(1:6), widths=chrom.lens$scaled, heights=rep(2, 6))
	par(mar=c(4, 4, 1, 1) + 0.1)
	
	for (chr in 1:6) {

		curr.win = win.smooth[[est]][[chr]]
		curr.roi = roi.smooth[[est]][roi.smooth[[est]]$chr == paste0("chr", chr),]
		
		# Get random window to show
		random.window.size = 100000000
		random.start = floor(runif(1, 1, 
								max(win.smooth[[est]][[chr]]$pos) - random.window.size))
		random.window = c(random.start, random.start + random.window.size)
		
		plot(curr.win$w.avg ~ curr.win$pos, type='n', 
			xlab=paste("Chromosome", chr), ylab=est.expr[[est]],
			ylim=ylims[[est]], xlim=random.window)

		# To plot the whole chromosome, use:
		# xlim=c(0, max(win.smooth[[est]][[1]]$pos))
		
		sig.cols = list()
		sig.cols[['50000']]  = 'red'
		sig.cols[['100000']] = 'purple'
		sig.cols[['150000']] = 'blue'
		
		for (sigma in sigmas) {
			this.sigma = curr.roi[curr.roi$sigma == sigma,]
			points(this.sigma$ests.wavg ~ this.sigma$middle, 
				type='p', pch='-', col=alpha(sig.cols[[sigma]], 0.5))
		}

		points(curr.win$w.avg ~ curr.win$pos, type='l', lty=1)
			
		# Add genes
		gene.height = ylims[[est]][1] + (0.8 * (ylims[[est]][2] - ylims[[est]][1]))
		points(x = genes[[chr]]$middle, y=rep(gene.height, nrow(genes[[chr]])), 
			cex=unlist(genes[[chr]]$gene.count) / 4,
			pch=16)
		
	}
	
	dev.off()
}

# ========================================================================================
# --- Output table of most extreme outliers
# ========================================================================================

bed = read.table("refGene/refGene.sort.gtf.papAnu2.bed", header=TRUE)

# ----------------------------------------------------------------------------------------
# --- Get most extreme outliers
# ----------------------------------------------------------------------------------------

for (est in ests) {

	extreme.hi = list()
	extreme.lo = list()
	for (sigma in sigmas) {
	
		curr.roi = roi.smooth[[est]][roi.smooth[[est]]$sigma == sigma,]
		curr.roi = merge(curr.roi, bed, by.x="start", by.y="Begin")
		curr.roi$Location = paste0(curr.roi$chr, ":", curr.roi$start, "-", curr.roi$end)
		
		rows.to.show = 20
		extreme.hi[[sigma]] = head(	curr.roi[order(curr.roi$p.val, decreasing=TRUE),], 
									n=rows.to.show)
		extreme.lo[[sigma]] = head(	curr.roi[order(curr.roi$p.val),],
									n=rows.to.show)
	}
	
	# Combine extremes of all sigmas
	extreme.hi = do.call(rbind, extreme.hi)
	extreme.lo = do.call(rbind, extreme.lo)
	
	# Set all but first instance of sigma value to ditto marks
	extreme.hi$sigma[-c(1 + 0:2 * rows.to.show)] = '"'
	extreme.lo$sigma[-c(1 + 0:2 * rows.to.show)] = '"'

	included.cols = c(10, 11, 6, 4, 5)
	display = c('s','s','s','s','f','f')
	digits  = c( 0,  0,  0,  0,  4,  4 )
	
	names(extreme.hi) = names(extreme.lo) = c("start", "chr", "end", "Weighted Estimate", 
												"$p$", "$\\sigma$", "middle", "Chr", 
												"End", "Gene", "Location")
	
	# Print high p-value table
	tex.out = paste0("results/ROI_smoothed_outliers_", 
							est, "_highp.tex")
	
	xt = xtable(extreme.hi[,included.cols],
					display=display, digits=digits)
	align(xt) = c('l', 'p{1in}||', rep('r', 4))

	sink(tex.out)
	print.xtable(xt, 
				include.rownames = FALSE, 
				sanitize.text.function = function(x){x},
				size="scriptsize",
				hline.after = c(-1, 0, 1:3 * rows.to.show))
	sink()

	# Print low p-value table
	tex.out = paste0("results/ROI_smoothed_outliers_", 
						est, "_lowp.tex")
	xt = xtable(extreme.lo[,included.cols],
					display=display, digits=digits)
	align(xt) = c('l', 'p{1in}||', rep('r', 4))

	sink(tex.out)
	print.xtable(xt, 
				include.rownames = FALSE, 
				sanitize.text.function = function(x){x},
				size="scriptsize",
				hline.after = c(-1, 0, 1:3 * rows.to.show))
	sink()
}

