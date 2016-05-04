#!/usr/bin/Rscript

# Call with:
# Rscript kernel_smooth_around_ROI.R [input.bed] [sigma] [estimates_to_smooth.txt]

library("snow")

cl = getMPIcluster()

options(stringsAsFactors = FALSE)
options(scipen=999)

# ------------------------------------------------------------------------------
# --- Get file of regions of interest
# ------------------------------------------------------------------------------

args = commandArgs(TRUE)

# Display usage if less than four arguments passed
if (length(args) < 4) {
  args = c("--help")
}

if ("--help" %in% args) {
	cat("
        Script for smoothing estimates, such as Fst, around regions of interest, such as genes

        Arguments:
        --help             - Show this text

		Usage:
			./kernel_smooth_around_ROI.R [roi.bed] [sigma] [estimates_to_smooth.txt] [label]
		
        Example to get Fst near testis-specific genes:
        	./kernel_smooth_around_ROI.R \
        		results/Macaque_Ensembl57_TopHat_UniqueReads.papAnu2.bed \
                        150000 \
                        results/all.pure.weir.fst \
                        fst\n\n")
 
		q(save="no")
}

# Get file of regions of interest, such as tissue-specific genes or all refGenes
roi.file = args[1]
roi = read.table(roi.file, header=TRUE)

### TEMPORARY FOR TESTING PURPOSES
###roi = head(roi, n=200)
###

sigma = as.numeric(args[2])		# e.g. 150000

est.file = as.character(args[3])	# File of estimates, e.g. Fst

out.label = as.character(args[4])	# Label for estimate, e.g. fst

# ------------------------------------------------------------------------------
# --- Get estimates to smooth (e.g. Fst)
# ------------------------------------------------------------------------------

ests = read.table(est.file, header=TRUE)

# Standardize header
names(ests) = c("chr", "pos", "estimate")

# Add chr prefix to chromosome if it is only a number
if (grepl("chr", ests[1]) == FALSE) {
   ests$chr = paste0("chr", ests$chr)
}

# ------------------------------------------------------------------------------
# --- Remove NaN values from estimates (Fst results, in particular)
# ------------------------------------------------------------------------------

ests = ests[ests$estimate != "NaN",]
ests = ests[is.na(ests$estimate) == FALSE,]

orig.estimates = list()
orig.estimates[[out.label]] = ests

# ------------------------------------------------------------------------------
# --- Function to get weighted average of stat given region
# --- and estimate p value via bootstrapping
# ------------------------------------------------------------------------------

compute.w.avg = function(chr, region.start, region.end, nreps.initial, nreps.max) {

	write(paste("Processing window ", chr, ": ", region.start, " to ", 
					region.end, "...", sep=""), stderr())

	# ------------------------------------------------------------------------------
	# --- Get weighted average of stat given region center
	# ------------------------------------------------------------------------------
	
	# Get all estimates in given range +/- 3*sigma
	current.ests = ests[ests$chr == chr &
						ests$pos >  (region.start - 3*sigma) &
						ests$pos <= (region.end   + 3*sigma),]
	
	# Remove any estimates that are NaN
	current.ests = current.ests[current.ests$estimate != "NaN",]
	
	# Abort if no SNPs are in this region
	if (dim(current.ests)[1] == 0) {
		write(paste("- No SNPs in window at location ", chr, ":", 
					region.start, "-", region.end, ". Skipping.", sep=""), stderr())
		return(c(chr, region.start, region.end, NA, NA))
	}
	
	# Weight by: exp[ -d^2 / (2*sigma^2) ]
	# where d is distance to nearest end of ROI
	
	weights.lo = exp( -(current.ests[current.ests$pos < region.start,]$pos 
			- region.start)^2 / (2*sigma^2) )
	weights.hi = exp( -(current.ests[current.ests$pos > region.end,]$pos 
			- region.end)^2 / (2*sigma^2) )
	n.snps.ongene = dim(current.ests)[1] - (length(c(weights.lo, weights.hi)))
	weights.ongene = rep(1, n.snps.ongene)

	weights = c(weights.lo, weights.ongene, weights.hi)
	
	#write(paste("position", chr, ":", current.ests$pos, 
	#			"has weight", weights), stderr())
	
	this.wavg = weighted.mean(current.ests$estimate, weights)
		
	#write(paste("Weighted avg of region centered at ", chr, ":", center.c, 
	#				"is", this.wavg), stderr())
	
	# ------------------------------------------------------------------------------
	# --- Estimate p value via bootstrapping
	# ------------------------------------------------------------------------------
	
	# Count of bootstrap averages that are greater than the current region's estimate
	gt.count = 0
	
	# Sample size
	sample.size = dim(current.ests)[1]
	
	#write(paste("Sampling with replacement. N =", sample.size), stderr())
	
	nreps = 0
	
	while (nreps <= nreps.max) {
			
		# If we have passed through the minimum number of reps
		if (nreps >= nreps.initial) {
					
			# And if the p-value is between 0.05 and 0.95 (boring)
			if ((((nreps - gt.count) / nreps) > 0.05) & (((nreps - gt.count) / nreps) < 0.95)) {
						
				# Quit with current p-value, since it seems uninteresting
				# and we have done the minimum number of replicates

				p.val = (nreps - gt.count) / nreps

				write(paste("Boring p-value @", nreps), stderr())	
				write(paste(" - p-value is", p.val), stderr())

				return(c(chr, region.start, region.end, this.wavg, p.val))
			
			}
		}
		
	
		# Sample the population genetic estimate, with replacement
		sample.wavg = weighted.mean(sample(	orig.estimates[[out.label]]$estimate, 
											sample.size, replace=TRUE), weights)
		
		#write(paste("Weighted avg of SAMPLE for region from", 
		#				chr, ":", region.start, "-", region.end, 
		#				"is", sample.wavg), stderr())
		
		if (this.wavg > sample.wavg) {
			gt.count = gt.count + 1
		}
		
		nreps = nreps + 1
	}
		
	p.val = (nreps - gt.count) / nreps
	
	#write (paste("Final p-value of", p.val, "after", nreps), stderr())

	return(c(chr, region.start, region.end, this.wavg, p.val))
}

# ------------------------------------------------------------------------------
# --- Compute weighted averages, p values for ROI
# ------------------------------------------------------------------------------

# Store results in big data.frame
# Chromosome, position.start, position.end, weighted average
all.weighted.avgs = list()

empty = data.frame(	chr = integer(0), start = character(0), end = character(0), 
					w.avg = character(0), p.val = character(0))

all.weighted.avgs[[out.label]] = empty

# Start with small number of replicants, step up until tails are accurate
#	      initial         max
#	Fst    10,000  10,000,000
	
# Compute weighted averages in range +/- 3*sigma for these estimates:

write(paste("- ", out.label, date()), stderr())
ptm <- proc.time()

#clusterExport(cl, varlist=c("compute.w.avg", "sigma", "roi", "ests", "orig.estimates"))
clusterExport(cl, c("compute.w.avg", "sigma", "roi", "ests", "orig.estimates", "out.label"))

###mc.cores = getOption("mc.cores", 2L)
###print(paste("Cores available:", mc.cores), stderr())

ests.wavg = parLapply(cl, 1:nrow(roi), function(x) {
	compute.w.avg(chr=roi[x,][["Chr"]], 
	region.start=roi[x,][["Begin"]], 
	region.end=roi[x,][["End"]], 
	nreps.initial=10000, nreps.max=10000000)
})

ests.wavg = as.data.frame(do.call(rbind, ests.wavg))
names(ests.wavg) = c("chr", "start", "end", "ests.wavg", "p.val")
write(proc.time() - ptm, stderr())
	
write(paste("- Finished ", date()), stderr())

# ------------------------------------------------------------------------------
# --- Output weighted average data files
# ------------------------------------------------------------------------------

out.suffix = paste(".", out.label, ".sigma", sigma, ".txt", sep="")
ests.out = gsub(".bed", out.suffix, roi.file)

write.table(ests.wavg,
	file=ests.out, 
	quote=FALSE, row.names=FALSE)

# ------------------------------------------------------------------------------
# --- Count significant values
# ------------------------------------------------------------------------------

low.pval.count  = nrow(ests.wavg[is.na(ests.wavg$p.val) == FALSE & ests.wavg$p.val < 0.05,])
high.pval.count = nrow(ests.wavg[is.na(ests.wavg$p.val) == FALSE & ests.wavg$p.val > 0.95,])

write(paste("Number of low  p-values (p < 0.05):", high.pval.count), stderr())
write(paste("Number of high p-values (p > 0.95):", low.pval.count),  stderr())

stopCluster(cl)
