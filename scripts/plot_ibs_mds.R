#!/usr/bin/Rscript

options(stringsAsFactors = FALSE)

# Generate individual_info.csv from the file on Google Drive
# Sequenced Awash Animals
# https://docs.google.com/spreadsheet/ccc?key=0AnJtKm_03RZDdEZkNEg4a2dTYkNjN1JvQlFUUFVhQ2c&usp=drive_web#gid=0

ind_info = read.csv("data/individual_info.csv")

mds = read.table('results/all.cleaned.LDpruned.mds', header = T)

group  = as.character(mds$FID)
taxon  = as.character(mds$FID)
ind_id = as.character(mds$FID)
csf    = as.character(mds$FID)
hyb_in = as.character(mds$FID)

for (i in 1:dim(ind_info)[1]) {

	this.ind = ind_info$Individual.ID[i]
	this.ind.info = ind_info[ind_info$Individual.ID == this.ind,]

	group  = replace(group,  group==this.ind, this.ind.info$Group)
	taxon  = replace(taxon,  taxon==this.ind, this.ind.info$Taxon)
	ind_id = replace(ind_id, ind_id==this.ind, this.ind.info$Sample.ID)
	csf    = replace(csf,    csf==this.ind, this.ind.info$CSF)
	hyb_in = replace(hyb_in, hyb_in==this.ind, this.ind.info$Hybrid.Index)

	# Eventually do same for Hybrid.Index, Sex, and Age.Estimate. Maybe Plate too?
}

# Replace "Filoha" with "X", so it doesn't become "F"
group[group == "Filoha"] = "X"

mds$group = group
mds$taxon = factor(taxon)
mds$ind_id = ind_id
mds$csf = csf
mds$hybrid_index = hyb_in

pdf(file='results/ibs_mds_plot.pdf')

	plot(mds$C1, mds$C2, pch = mds$group, col=mds$taxon, 
		main="IBS MDS Plot",
		xlab="Dimension 1", ylab="Dimension 2", asp=1)

	legend("bottomright", inset=.05,
		c("Anubis","Hamadryas", "Hybrids"), 
		fill=1:4, horiz=FALSE)

dev.off()

# Export mds now that info on individuals has been added.
write.table(mds, file="results/mds.txt")
