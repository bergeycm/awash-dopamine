# ------------------------------------------------------------------------------
# --- Plot moving averages. Find and shade outliers
# ------------------------------------------------------------------------------

library(ggplot2)

options(stringsAsFactors = TRUE)
options(scipen=999)

#alpha = 0.001	# p < 10^(-3)
sigma       = 150000
chr.to.show = 1:20

# Bring in info on spermatogenesis genes to add to plot
sperm = read.table("data/spermatogenesis_genes.txt", header=TRUE)
# Sort it numerically
sperm$chr = gsub("chr", "", sperm$chr)
sperm = sperm[order(as.numeric(sperm$chr), sperm$start),]
sperm$chr = paste('chr', sperm$chr, sep='')

# Bring in info on MHC genes to add to plot
mhc = read.table("data/hla_baboon_orthologs.txt", header=TRUE)
# Sort it numerically
mhc$chr = gsub("chr", "", mhc$chr)
mhc = mhc[order(as.numeric(mhc$chr), mhc$start),]
mhc$chr = paste('chr', mhc$chr, sep='')

# Read in data output from kernel_smooth.R
all.weighted.avgs = list()

empty = data.frame(chr = integer(0), pos = character(0), w.avg = character(0), p.val = character(0))

# Skipping all but Fst for now
all.weighted.avgs[["PI_HAMA"]]  = empty
all.weighted.avgs[["PI_ANUB"]]  = empty
all.weighted.avgs[["HET_HAMA"]] = empty
all.weighted.avgs[["HET_ANUB"]] = empty
all.weighted.avgs[["FST"]]      = empty

for (chr in paste("chr", chr.to.show, sep="")) {
	all.weighted.avgs[["FST"]] = rbind(all.weighted.avgs[["FST"]], 
									read.table(paste("results/fst_weighted_averages_", 
									chr, ".txt", sep=""), header=TRUE))
}

plot.along.chrom = function(est.name, pdf.prefix, y.label, alpha=0.01, plot.high=TRUE, plot.low=TRUE, add.sperm=FALSE) {

	outliers.lo = all.weighted.avgs[[est.name]][all.weighted.avgs[[est.name]]$p.val > (1-alpha),]
	outliers.hi = all.weighted.avgs[[est.name]][all.weighted.avgs[[est.name]]$p.val < alpha,]

	outliers.lo.coords = data.frame(chr=outliers.lo$chr, start=outliers.lo$pos - 2 * sigma, end=outliers.lo$pos + 2 * sigma)
	outliers.hi.coords = data.frame(chr=outliers.hi$chr, start=outliers.hi$pos - 2 * sigma, end=outliers.hi$pos + 2 * sigma)
	
	this.est = all.weighted.avgs[[est.name]]

	q = ggplot(this.est)

	# Add ymax and ymin to these outlier coordinates
	# Hacky way to get ggplot to access it
	outliers.lo.coords$this.ymin = min(this.est$w.avg, na.rm=TRUE)
	outliers.lo.coords$this.ymax = max(this.est$w.avg, na.rm=TRUE)
	outliers.hi.coords$this.ymin = min(this.est$w.avg, na.rm=TRUE)
	outliers.hi.coords$this.ymax = max(this.est$w.avg, na.rm=TRUE)

	# Add low outliers
	if (plot.low) {
		q = q + geom_rect(	data=outliers.lo.coords, 
							mapping=aes(xmin=start,     xmax=end, 
										ymin=this.ymin, ymax=this.ymax), 
							fill="yellow", alpha=0.2)
	}
	
	# Add high outliers
	if (plot.high) {
		q = q + geom_rect(	data=outliers.hi.coords, 
							mapping=aes(xmin=start,     xmax=end, 
										ymin=this.ymin, ymax=this.ymax), 
							fill="orange", alpha=0.2)
	}
	
	# Add spermatogenesis genes
	if (add.sperm) {
		q = q + geom_vline(aes(xintercept = (sperm$start + sperm$end) / 2), sperm, col="red")
	}
	
	# Add MHC genes
	#q = q + geom_vline(aes(xintercept = (mhc$start + mhc$end) / 2), mhc, col="blue", alpha=0.2)
	
	q = q + geom_line(aes(x = pos, y = w.avg)) + theme_bw() +
				ylab(y.label) + facet_grid(chr ~ .)

	ggsave(filename=paste(pdf.prefix, ".pdf", sep=''), plot=q, height=12.5, width=25)
	
	write(paste("Saving ", est.name, " plot to: ", pdf.prefix, ".pdf", sep=''), stderr())

}

# Skipping all but Fst for now
#plot.along.chrom(est.name="PI_HAMA",  pdf.prefix="results/pi_hama_",  y.label=expression("Pure Hamadryas - " * pi))
#plot.along.chrom(est.name="PI_ANUB",  pdf.prefix="results/pi_anub_",  y.label=expression("Pure Anubis - " * pi))
#plot.along.chrom(est.name="HET_HAMA", pdf.prefix="results/het_hama_", y.label=expression("Pure Hamadryas - H"["O"]))
#plot.along.chrom(est.name="HET_ANUB", pdf.prefix="results/het_anub_", y.label=expression("Pure Anubis - H"["O"]))

# Plot extreme values of Fst
plot.along.chrom(	est.name="FST",
					pdf.prefix=paste("results/fst_0.001", sep=""),
					y.label=expression('F'["ST"]),
					alpha=0.001,
					plot.high=TRUE,
					plot.low=TRUE,
					add.sperm=FALSE)

# Show all regions used in PANTHER
plot.along.chrom(	est.name="FST",
					pdf.prefix=paste("results/fst_0.05", sep=""),
					y.label=expression('F'["ST"]),
					alpha=0.05,
					plot.high=TRUE,
					plot.low=TRUE,
					add.sperm=FALSE)

# Show high regions and spermatogenesis gene sites
plot.along.chrom(	est.name="FST",
					pdf.prefix=paste("results/fst_0.05_sperm", sep=""),
					y.label=expression('F'["ST"]),
					alpha=0.05,
					plot.high=TRUE,
					plot.low=FALSE,
					add.sperm=TRUE)
