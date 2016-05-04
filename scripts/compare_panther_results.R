#!/usr/bin/env Rscript

sigmas = c('50000', '100000', '150000')
ests   = c('fst')
annos  = c('go.bp', 'go.cc', 'go.mf', 'panther.proteinclass', 'pathway')

# Create directory for output files
out.dir = "results/panther_results/methods_comparison/"
dir.create(out.dir, showWarnings = FALSE)

# Function to plot p-values of different tests
plot.p.vals = function (dir=c('hi', 'lo')) {

	par(mfrow=c(total.cols,total.cols))
	par(pty="s")
	par(mar=c(rep(1, 4)))
	par(oma=c(rep(4, 4)))
	for (i in 1:total.cols) {
		for (j in 1:total.cols) {
		
			if (i == j) {
				# Print title
				plot(c(0, 1), c(0, 1), ann=FALSE, bty='n', type='n', 
					xaxt='n', yaxt='n')
				text(x = 0.5, y = 0.5,
					names(enrich.p.vals.hi)[i], 
					cex = 1.6, col = "black")
					
			} else {
				if (dir == 'hi') {
					these.vals = enrich.p.vals.hi
				} else if (dir == 'lo') {
					these.vals = enrich.p.vals.lo
				}
				plot(these.vals[,i], these.vals[,j], 
					main='', xlab='', ylab='', asp=1, axes=FALSE,
					xlim=c(0,1), ylim=c(0,1), xaxs='i', yaxs='i')
				box()
				
				# Add lines denoting significant p-values
				abline(h=0.05, col='red', lty=2)
				abline(v=0.05, col='red', lty=2)
			}

			# Add axis for only outside edges of graphs
			if (i == 1) {
				# Axes on top
				axis(side=3, labels=TRUE, tck=-0.02)
			} else if (i == total.cols) {
				# Axes on bottom
				axis(side=1, labels=TRUE, tck=-0.02)
			}
			
			if (j == 1) {
				# Axes on left
				axis(side=2, labels=TRUE, tck=-0.02)
			} else if (j == total.cols) {
				# Axes on right
				axis(side=4, labels=TRUE, tck=-0.02)
			}
		}
	}
	if (dir == 'hi') { sign = '+' }
	if (dir == 'lo') { sign = '-' }
	mtext(paste0(est, " - ", anno, " (", sign, ")"), side=3, line=2, outer=TRUE)
}

# For each estimate
for (est in ests) {

	# For each annotation type
	for (anno in annos) {

		# --------------------------------------------------------------------------------
		# --- Read in PANTHER output files
		# --------------------------------------------------------------------------------
				
		e_roi = lapply(	sigmas, 
						function (sigma) {
							read.table(
								paste0("refGene/refGene.sort.gtf.papAnu2.", est, 
										".sigma", sigma, ".panther.", est, 
										".panther.enrich.", anno, ".txt")
							)
						}
					)
		names(e_roi) = sigmas

		e_win = read.table(paste0(	"results/panther_results/", est, 
									".all.orthologs.genes.list.panther.enrich.", 
									anno, ".txt"))
		o_win = read.table(paste0(	"results/panther_results/", est, 
									".overrep.", anno, ".txt"), header=TRUE)

		# --------------------------------------------------------------------------------
		# --- Make data.frame of p-values
		# --------------------------------------------------------------------------------
		
		enrich.p.vals = as.data.frame(cbind(sapply(e_roi, function(x) x$p), e_win$p))
		names(enrich.p.vals) = c(paste("ROI", sigmas), "window")

		# Make copy reduced to just overrepresented or underrepresented annotations
		enrich.p.vals.hi = enrich.p.vals
		enrich.p.vals.lo = enrich.p.vals

		# Switch ROI values to NA if they do not have correct sign
		total.cols = ncol(enrich.p.vals.hi)
		for (i in 1:(total.cols - 1)) {
			enrich.p.vals.hi[[i]] [e_roi[[i]]$overUnder != "+"] = NA
			enrich.p.vals.lo[[i]] [e_roi[[i]]$overUnder != "-"] = NA
		}
		# Switch windowed values to NA if they do not have correct sign
		enrich.p.vals.hi[[total.cols]] [e_win$overUnder != "+"] = NA
		enrich.p.vals.lo[[total.cols]] [e_win$overUnder != "-"] = NA

		output.prefix = paste0(out.dir, "methods_comparison_", est, "_", anno)
		pdf(paste0(output.prefix, ".pdf"))
			plot.p.vals('hi')
			plot.p.vals('lo')
		dev.off()
		
		enrich.p.vals.hi.noNA = enrich.p.vals.hi[complete.cases(enrich.p.vals.hi),]
		enrich.p.vals.lo.noNA = enrich.p.vals.lo[complete.cases(enrich.p.vals.lo),]
		
		sink(paste0(output.prefix, ".txt"))
			print(cor(enrich.p.vals.hi.noNA))
			print(cor(enrich.p.vals.lo.noNA))
		sink()
	}
}
