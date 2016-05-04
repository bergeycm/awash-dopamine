#!/usr/bin/Rscript

# Call with:
# Rscript kernel_smooth.R chr1 --fst           # Just chr1, fst
# Rscript kernel_smooth.R chr18 --pi_a --fst   # Just chr18, pi hamadryas and fst
# Rscript kernel_smooth.R --fst                # Just fst, all chromosomes

library("parallel")
library("foreach")

cl = makeCluster(detectCores())

# ------------------------------------------------------------------------------
# --- Smooth estimates of nucleotide diversity, heterozygosity, and Fst
# ------------------------------------------------------------------------------

options(stringsAsFactors = FALSE)
options(scipen=999)

pi.hama = read.table("results/pure.hamadryas.sites.pi", header=TRUE)
pi.anub = read.table("results/pure.anubis.sites.pi",    header=TRUE)

hwe.hama = read.table("results/pure.hamadryas.hwe", header=TRUE)
hwe.anub = read.table("results/pure.anubis.hwe",    header=TRUE)

fst = read.table("results/all.pure.weir.fst", header=TRUE)

taj.d.hama = read.table("results/pure.hamadryas.Tajima.D", header=TRUE)
taj.d.anub = read.table("results/pure.anubis.Tajima.D",    header=TRUE)

sigma     = 150000
step.size = 100000

# ------------------------------------------------------------------------------
# --- Process arguments from user (chromosome and analyses requested)
# ------------------------------------------------------------------------------

args = commandArgs(TRUE)

# Display usage if no arguments passed
if (length(args) < 1) {
  args = c("--help")
}
 
if ("--help" %in% args) {
	cat("
        Kernel-smoothing, significance-testing script for population genetic statistics.

        Arguments:
        --pi_h             - Do pi for hamadryas
        --pi_a             - Do pi for anubis
        --het_h            - Do heterozygosity for hamadryas
        --het_a            - Do heterozygosity for anubis
        --fst              - Do Fst
       --help             - Show this text

        Example to do pi for hamadryas and Fst:
        ./kernel_smooth.R chr18 --pi_h --fst\n\n")
 
		q(save="no")
}

do.pi.hama  = do.pi.anub  = FALSE
do.het.hama = do.het.anub = FALSE
do.fst = FALSE

if ("--pi_h" %in% args) {
	do.pi.hama = TRUE
}
if ("--pi_a" %in% args) {
	do.pi.anub = TRUE
}
if ("--het_h" %in% args) {
	do.het.hama = TRUE
}
if ("--het_a" %in% args) {
	do.het.anub = TRUE
}
if ("--fst" %in% args) {
	do.fst = TRUE
}

# ------------------------------------------------------------------------------
# --- Get chromosomes to process from user. Default is to do all.
# ------------------------------------------------------------------------------

# Anything that matches "chr\\d+" will be considered a chromosome to do
chrs.to.do = grep("chr\\d+", args, perl=TRUE, value=TRUE)

if (length(chrs.to.do) == 0) {
	chrs.to.do = row.names(table(pi.hama$CHROM))
}

write(paste("Processing chromosome(s):", chrs.to.do), stderr())

# ------------------------------------------------------------------------------
# --- Remove zero values from pi estimates
# ------------------------------------------------------------------------------

pi.hama = pi.hama[pi.hama$PI != 0,]
pi.anub = pi.anub[pi.anub$PI != 0,]

# ------------------------------------------------------------------------------
# --- Remove NaN values from fst
# ------------------------------------------------------------------------------

fst = fst[fst$WEIR_AND_COCKERHAM_FST != "NaN",]

# ------------------------------------------------------------------------------
# --- Fix columns of HWE data.frame
# ------------------------------------------------------------------------------

# Split the columns of hwe.* that have headers like:
# 	OBS.HOM1.HET.HOM2.
# and values like:
# 	15/6/2

names(hwe.hama) = names(hwe.anub) = c("CHROM", "POS", "OBS", "EXP", "ChiSq", "P")

hwe.hama = within(hwe.hama, OBS<-data.frame(do.call('rbind', strsplit(as.character(OBS), '/', fixed=TRUE))))
hwe.hama = within(hwe.hama, EXP<-data.frame(do.call('rbind', strsplit(as.character(EXP), '/', fixed=TRUE))))

hwe.anub = within(hwe.anub, OBS<-data.frame(do.call('rbind', strsplit(as.character(OBS), '/', fixed=TRUE))))
hwe.anub = within(hwe.anub, EXP<-data.frame(do.call('rbind', strsplit(as.character(EXP), '/', fixed=TRUE))))

names(hwe.hama$OBS) = names(hwe.hama$EXP) = names(hwe.anub$OBS) = names(hwe.anub$EXP) = c("HOM1", "HET", "HOM2")

# Fix class
hwe.hama.obs = as.data.frame(lapply(hwe.hama[,"OBS"], as.numeric))
hwe.anub.obs = as.data.frame(lapply(hwe.anub[,"OBS"], as.numeric))

# Compute percentage of animals that are heterozygosites (hwe.hama[,"OBS"][2]).
# Put value to its own column and convert it to numeric class
hwe.hama[,"HET"] = as.numeric(unlist(hwe.hama[,"OBS"][2])) / rowSums(hwe.hama.obs)
hwe.anub[,"HET"] = as.numeric(unlist(hwe.anub[,"OBS"][2])) / rowSums(hwe.anub.obs)

# Get rid of zero values of heterozygosity
hwe.hama = hwe.hama[hwe.hama[,"HET"] != 0,]
hwe.anub = hwe.anub[hwe.anub[,"HET"] != 0,]

# ------------------------------------------------------------------------------
# --- Make list to store original (unsmoothed) values
# ------------------------------------------------------------------------------

orig.estimates = list()
orig.estimates[["pi.hama"]] = pi.hama
orig.estimates[["pi.anub"]] = pi.anub
orig.estimates[["hwe.hama"]] = hwe.hama
orig.estimates[["hwe.anub"]] = hwe.anub
orig.estimates[["fst"]] = fst

# ------------------------------------------------------------------------------
# --- Function to get weighted average of stat given region center
# --- and estimate p value via bootstrapping
# ------------------------------------------------------------------------------

compute.w.avg = function(pop.gen.ests, col.name, orig.est.df, center.c, nreps.initial, nreps.max) {

	write(paste("Processing window", center.c, "...", sep=""), stderr())

	# ------------------------------------------------------------------------------
	# --- Get weighted average of stat given region center
	# ------------------------------------------------------------------------------
	
	# Get all estimates in range center.c +/- 3*sigma
	current.ests = pop.gen.ests[pop.gen.ests$CHROM == chr &
								pop.gen.ests$POS >  (center.c - 3*sigma) &
								pop.gen.ests$POS <= (center.c + 3*sigma),]
	
	# Remove any estimates that are NaN
	current.ests = current.ests[current.ests[,col.name] != "NaN",]
	
	# Abort if no SNPs are in this region
	if (dim(current.ests)[1] == 0) {
		write(paste("- No SNPs in window at location", chr, ":", current.ests$POS, 
					". Skipping.", sep=""), stderr())
		return(c(center.c, NA, NA))
	}
	
	# Weight by: exp[ -(p-c)^2 / (2*sigma^2) ]
	weights = exp( -(current.ests$POS - center.c)^2 / (2*sigma^2) )

	#write(paste("position", chr, ":", current.ests$POS, 
	#			"has weight", weights), stderr())
	
	this.wavg = weighted.mean(current.ests[,col.name], weights)
		
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

				#write(paste("Boring p-value @", nreps), stderr())	
				#write(paste(" - p-value is", p.val), stderr())

				return(c(center.c, this.wavg, p.val))
			
			}
		}
		
	
		# Sample the population genetic estimate, with replacement
		sample.wavg = weighted.mean(sample(orig.estimates[[orig.est.df]][,col.name], sample.size, replace=TRUE), weights)
		
		#write(paste("Weighted avg of SAMPLE for region centered at ", chr, ":", center.c, 
		#				"is", sample.wavg), stderr())
		
		if (this.wavg > sample.wavg) {
			gt.count = gt.count + 1
		}
		
		nreps = nreps + 1
	}
		
	p.val = (nreps - gt.count) / nreps
	
	#write (paste("Final p-value of", p.val, "after", nreps), stderr())

	return(c(center.c, this.wavg, p.val))
}

# ------------------------------------------------------------------------------
# --- Compute weighted averages, p values along each chromosome
# ------------------------------------------------------------------------------

# Store results in big data.frame
# Chromosome, position, weighted average
all.weighted.avgs = list()

empty = data.frame(chr = integer(0), pos = character(0), w.avg = character(0), p.val = character(0))

all.weighted.avgs[["PI_HAMA"]]  = empty
all.weighted.avgs[["PI_ANUB"]]  = empty
all.weighted.avgs[["HET_HAMA"]] = empty
all.weighted.avgs[["HET_ANUB"]] = empty
all.weighted.avgs[["FST"]]      = empty

for (chr in chrs.to.do) {

	write(paste("Analyzing chromosome:", chr), stderr())
	
	# Get highest SNP position on this chromosome
	max.snp = max(pi.hama$POS[pi.hama$CHROM == chr])

	# Go from 1 to the approx. length of the chromosome, in steps of size step.size
	window.starts = seq(from=1, to=max.snp, by=step.size)
	window.middles = window.starts + (0.5 * step.size)
	
	# Start with small number of replicants, step up until tails are accurate
	#	          initial         max
	#	pi,H          100   1,000,000
	#	D           1,000   1,000,000
	#	Fst        10,000  10,000,000
	
	# Compute weighted averages in range +/- 3*sigma
	# ...for pi:
	if (do.pi.hama) {
		write(paste("- Pi - Hamadryas", date()), stderr())
		ptm <- proc.time()
		pi.hama.wavg = do.call(rbind, lapply(window.middles, compute.w.avg, 
			pop.gen.ests=pi.hama, col.name="PI", orig.est.df="pi.hama", 
			nreps.initial=100, nreps.max=1000000))
		pi.hama.wavg = as.data.frame(pi.hama.wavg)
		names(pi.hama.wavg) = c("pos", "pi.wavg", "p.val")
		write(proc.time() - ptm, stderr())
	}
	
	if (do.pi.anub) {
		write(paste("- Pi - Anubis", date()), stderr())
		ptm <- proc.time()
		pi.anub.wavg = do.call(rbind, lapply(window.middles, compute.w.avg, 
			pop.gen.ests=pi.anub, col.name="PI", orig.est.df="pi.anub", 
			nreps.initial=100, nreps.max=1000000))
		pi.anub.wavg = as.data.frame(pi.anub.wavg)
		names(pi.anub.wavg) = c("pos", "pi.wavg", "p.val")
		write(proc.time() - ptm, stderr())
	}
	
	# ...for observed heterozygosity:
	if (do.het.hama) {
		write(paste("- Heterozygosity - Hamadryas", date()), stderr())
		ptm <- proc.time()
		hwe.hama.wavg = do.call(rbind, lapply(window.middles, compute.w.avg, 
			pop.gen.ests=hwe.hama, col.name="HET", orig.est.df="hwe.hama", 
			nreps.initial=100, nreps.max=1000000))
		hwe.hama.wavg = as.data.frame(hwe.hama.wavg)
		names(hwe.hama.wavg) = c("pos", "hwe.wavg", "p.val")
		write(proc.time() - ptm, stderr())
	}
	
	if (do.het.anub) {
		write(paste("- Heterozygosity - Anubis", date()), stderr())
		ptm <- proc.time()
		hwe.anub.wavg = do.call(rbind, lapply(window.middles, compute.w.avg, 
			pop.gen.ests=hwe.anub, col.name="HET", orig.est.df="hwe.anub", 
			nreps.initial=100, nreps.max=1000000))
		hwe.anub.wavg = as.data.frame(hwe.anub.wavg)
		names(hwe.anub.wavg) = c("pos", "hwe.wavg", "p.val")
		write(proc.time() - ptm, stderr())
	}	
	
	# ...for Fst:
	if (do.fst) {
		write(paste("- Fst", date()), stderr())
		ptm <- proc.time()
		clusterExport(cl, varlist=c("chr", "sigma", "orig.estimates"))
		fst.wavg = do.call(rbind, parLapply(cl, window.middles, compute.w.avg, 
			pop.gen.ests=fst, col.name="WEIR_AND_COCKERHAM_FST", orig.est.df="fst", 
			nreps.initial=10000, nreps.max=10000000))
		fst.wavg = as.data.frame(fst.wavg)
		names(fst.wavg) = c("pos", "fst.wavg", "p.val")
		write(proc.time() - ptm, stderr())
	}
	
	write(paste("- Finished ", date()), stderr())

	# ------------------------------------------------------------------------------
	# --- Add this chromosome's measurements to big data.frame
	# ------------------------------------------------------------------------------
	
	if (do.pi.hama) {
		all.weighted.avgs$PI_HAMA = rbind(all.weighted.avgs$PI_HAMA, 
				setNames(cbind(chr, pi.hama.wavg), names(empty)))
	}
	if (do.pi.anub) {
		all.weighted.avgs$PI_ANUB = rbind(all.weighted.avgs$PI_ANUB, 
				setNames(cbind(chr, pi.anub.wavg), names(empty)))
	}
	if (do.het.hama) {
		all.weighted.avgs$HET_HAMA = rbind(all.weighted.avgs$HET_HAMA, 
				setNames(cbind(chr, hwe.hama.wavg), names(empty)))
	}
	if (do.het.anub) {
		all.weighted.avgs$HET_ANUB = rbind(all.weighted.avgs$HET_ANUB, 
				setNames(cbind(chr, hwe.anub.wavg), names(empty)))
	}
	if (do.fst) {
		all.weighted.avgs$FST = rbind(all.weighted.avgs$FST, 
				setNames(cbind(chr, fst.wavg), names(empty)))
	}
	
}

# ------------------------------------------------------------------------------
# --- Find most extreme Tajima's D values
# ------------------------------------------------------------------------------

# Alpha value for extreme values of Tajima's D
tajd.alpha = 0.001
tajd.probs = c(tajd.alpha, 1 - tajd.alpha)

# Find cut-off values
tajd.cutoffs.hama = quantile(taj.d.hama[taj.d.hama$N_SNPS > 0,]$TajimaD, 
								probs=tajd.probs)
tajd.cutoffs.anub = quantile(taj.d.anub[taj.d.anub$N_SNPS > 0,]$TajimaD, 
								probs=tajd.probs)

tajd.outliers.lo.hama = taj.d.hama[taj.d.hama$TajimaD < tajd.cutoffs.hama[1],]
tajd.outliers.hi.hama = taj.d.hama[taj.d.hama$TajimaD > tajd.cutoffs.hama[2],]

tajd.outliers.lo.anub = taj.d.anub[taj.d.anub$TajimaD < tajd.cutoffs.anub[1],]
tajd.outliers.hi.anub = taj.d.anub[taj.d.anub$TajimaD > tajd.cutoffs.anub[2],]


# ------------------------------------------------------------------------------
# --- Output weighted average data files
# ------------------------------------------------------------------------------

all.chrs = paste(chrs.to.do, sep="-", collapse="-")

if (do.pi.hama) {
	write.table(all.weighted.avgs$PI_HAMA,
		file=paste("results/pi-hama_weighted_averages_",  all.chrs, ".txt", sep=""), 
		quote=FALSE, row.names=FALSE)
}
if (do.pi.anub) {
	write.table(all.weighted.avgs$PI_ANUB,
		file=paste("results/pi-anub_weighted_averages_",  all.chrs, ".txt", sep=""), 
		quote=FALSE, row.names=FALSE)
}
if (do.het.hama) {
	write.table(all.weighted.avgs$HET_HAMA,
		file=paste("results/het-hama_weighted_averages_", all.chrs, ".txt", sep=""), 
		quote=FALSE, row.names=FALSE)
}
if (do.het.anub) {
	write.table(all.weighted.avgs$HET_ANUB,
		file=paste("results/het-anub_weighted_averages_", all.chrs, ".txt", sep=""), 
		quote=FALSE, row.names=FALSE)
}
if (do.fst) {
	write.table(all.weighted.avgs$FST,
		file=paste("results/fst_weighted_averages_",      all.chrs, ".txt", sep=""), 
		quote=FALSE, row.names=FALSE)
}

stopCluster(cl)
