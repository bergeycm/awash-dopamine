#!/usr/bin/Rscript

options(stringsAsFactors = FALSE)

mds = read.table('results/mds.txt', header = T)

adm = read.table("results/all.cleaned.LDpruned.2.Q")

# ========================================================================================
# === Make lists of pure anubis and hamadryas animals
# ========================================================================================

#group = mds$group
#ind_id = mds$ind_id

#new.order = order(adm$V1, group, ind_id)

# Give the ADMIXTURE results descriptive names
names(adm) = c("ADMIXTURE_1", "ADMIXTURE_2")

# Add the ADMIXTURE results to the MDS results (with extra info)
mds.adm = cbind(mds, adm)

# Find extremes
pure_high_adm = mds.adm[mds.adm$ADMIXTURE_1 < 0.001,]
pure_low_adm  = mds.adm[mds.adm$ADMIXTURE_1 > 0.999,]

# Find hybrids
hybrids_adm = mds.adm[(mds.adm$ADMIXTURE_1 >= 0.001) & (mds.adm$ADMIXTURE_1 <= 0.999),]

# Figure out which represents the pure anubis and which the pure hamadyras
# since the ADMIXTURE results could be inverted

num.high.hama = length(pure_high_adm$taxon[pure_high_adm$taxon == "Hamadryas"])
num.high.anub = length(pure_high_adm$taxon[pure_high_adm$taxon == "Anubis"])

if (num.high.hama > num.high.anub) {
	# Subset with high ADMIXTURE values is predominantly Hamadryas
	pure.hama.subset = pure_high_adm
	pure.anub.subset = pure_low_adm
} else {
	# Subset with high ADMIXTURE values is predominantly Anubis
	pure.hama.subset = pure_low_adm
	pure.anub.subset = pure_high_adm
}

# Save these datasets as text files
write.table(pure.hama.subset, file="results/pure_hama_data.txt")
write.table(pure.anub.subset, file="results/pure_anub_data.txt")
write.table(hybrids_adm,      file="results/hybrid_data.txt")

# Also output first two columns (FID and IID, which are currently the same thing)
# for plink to read.

write.table(pure.hama.subset[,1:2], file="results/pure_hama_list.txt", 
	quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
write.table(pure.anub.subset[,1:2], file="results/pure_anub_list.txt", 
	quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
write.table(hybrids_adm[,1:2],      file="results/hybrid_list.txt", 
	quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
