#!/usr/bin/Rscript

options(stringsAsFactors = FALSE)

args = commandArgs(trailingOnly = TRUE)
prefix = args[1]
if (is.na(prefix)) {
	prefix=""
}

admix_cv_error = read.csv(paste("reports/", prefix, "ADMIXTURE_log_cv.out", sep=""), sep=":", header=FALSE)

pdf(file=paste('reports/', prefix, 'ADMIXTURE_CV_plot.pdf', sep=""))

	plot(admix_cv_error$V2 ~ seq(nrow(admix_cv_error)), type="b", pch=19,
		main="ADMIXTURE Cross-Validation Error Plot",
		xlab="K", ylab="Cross-Validation Error")
	
	ideal.k = 2
	points(admix_cv_error$V2[ideal.k] ~ ideal.k, col="red", cex=2)
	text(ideal.k, admix_cv_error$V2[ideal.k], 
		labels=round(admix_cv_error$V2[ideal.k], digits=3), 
		cex=0.8, pos=2, col="red")

	lowest.k = which (admix_cv_error$V2 == min(admix_cv_error$V2))
	points(admix_cv_error$V2[lowest.k] ~ lowest.k, col="red", cex=2)
	text(lowest.k, admix_cv_error$V2[lowest.k], 
		labels=round(admix_cv_error$V2[lowest.k], digits=3), 
		cex=0.8, pos=3, col="red")
		
dev.off()
