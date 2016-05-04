#!/usr/bin/Rscript

library (ggplot2)
library(car)

options(stringsAsFactors = FALSE)

mds = read.table('results/mds.txt', header = T)

# ========================================================================================
# === Analyze and plot ADMIXTURE results
# ========================================================================================

adm.files = list.files(path = "results/", pattern = "all.*.Q$")
pattern = "[0-9]+.Q"
k.vals = gsub(".Q", "", regmatches(adm.files, regexpr(pattern, adm.files)))

adm = list()
se  = list()

# Hard code colors
cols = list()
cols[["2"]] = c("#FF7A00", "#03899C")
cols[["3"]] = c("#FF7A00", "#A001A6", "#03899C")
cols[["6"]] = c("#FF7A00", "#58E000", "#FFC300", "#E8003A", "#2219B2", "#03899C")

for (k in k.vals) {

	k.ch = as.character(k)
	adm[[k.ch]] = read.table(paste0("results/all.cleaned.LDpruned.", k, ".Q"))

	# Swap in consistent order, based on total amount of ancestry ascribed to component
	col.order = order(colSums(adm[[k.ch]]))
	adm[[k.ch]] = adm[[k.ch]][,col.order]
	names(adm[[k.ch]]) = paste0("V", 1:ncol(adm[[k.ch]]))
	
	group = mds$group
	ind_id = mds$ind_id
	
	num.ind = nrow(mds)
	
	new.order = order(adm[[k.ch]]$V1,
		if ("V2" %in% names(adm[[k.ch]])) { adm[[k.ch]]$V2 } else { rep(NA, num.ind) },
		if ("V3" %in% names(adm[[k.ch]])) { adm[[k.ch]]$V3 } else { rep(NA, num.ind) },
		if ("V4" %in% names(adm[[k.ch]])) { adm[[k.ch]]$V4 } else { rep(NA, num.ind) },
		if ("V5" %in% names(adm[[k.ch]])) { adm[[k.ch]]$V5 } else { rep(NA, num.ind) },
		if ("V6" %in% names(adm[[k.ch]])) { adm[[k.ch]]$V6 } else { rep(NA, num.ind) },
		mds$group,
		mds$ind_id)
	
	# Plot animals by ancestry, largely ignoring group
	if (k == 2) {
		pdf(file='results/admixture_plot.pdf', width=15, height=5)
	
			barplot(t(as.matrix(adm[[k.ch]][new.order,])), 
				col=cols[[k.ch]], 
				xlab="Individual", ylab="Ancestry", border=NA, xaxt='n')
		
		dev.off()
	}
	
	# Plot animals by group first, then ancestry within
	new.order.group = order(mds$group,
		adm[[k.ch]]$V1,
		if ("V2" %in% names(adm[[k.ch]])) { adm[[k.ch]]$V2 } else { rep(NA, num.ind) },
		if ("V3" %in% names(adm[[k.ch]])) { adm[[k.ch]]$V3 } else { rep(NA, num.ind) },
		if ("V4" %in% names(adm[[k.ch]])) { adm[[k.ch]]$V4 } else { rep(NA, num.ind) },
		if ("V5" %in% names(adm[[k.ch]])) { adm[[k.ch]]$V5 } else { rep(NA, num.ind) },
		if ("V6" %in% names(adm[[k.ch]])) { adm[[k.ch]]$V6 } else { rep(NA, num.ind) },
		mds$ind_id)

	groups = mds[new.order.group,]$group
	group.starts = which(groups != c(groups[-1], NA))

	# Read in file with SE estimates
	se[[k.ch]] = read.table(paste0("results/all.cleaned.LDpruned.", k, ".Q_se"))
	
	pdf(file=paste0('results/admixture_plot_K', k, '_by_group.pdf'), width=15, height=5)
	
		mp = barplot(t(as.matrix(adm[[k.ch]][new.order.group,])), 
			col=cols[[k.ch]], 
			xlab="Group", ylab="Ancestry", border=NA, xaxt='n')
	
		bar.width = mp[2] - mp[1]
		border.locs = mp[group.starts] + (0.5 * bar.width)
		abline(v=border.locs, lwd=3)

		group.mps = colMeans(rbind(	c(0, border.locs), 
									c(border.locs, max(mp) + 0.5*bar.width)))

		uniq.groups = unique(groups)
	
		mtext(uniq.groups, side=1, at=group.mps, padj=1)

		if (k == 2) {
			arrows(	mp, adm[[k.ch]][new.order.group,]$V1+se[[k.ch]][new.order.group,]$V1, 
					mp, adm[[k.ch]][new.order.group,]$V1, 
					angle=90, code=1, length=0)
			arrows(	mp, adm[[k.ch]][new.order.group,]$V1-se[[k.ch]][new.order.group,]$V1, 
					mp, adm[[k.ch]][new.order.group,]$V1, 
					angle=90, code=1, length=0)
		}
	dev.off()

}

# ========================================================================================
# === Plot various ancestry inferences against one another
# ========================================================================================

# Set adm equal to k=2 run, just for simplicity's sake
adm = read.table(paste0("results/all.cleaned.LDpruned.2.Q"))

# Swap ADM1 and ADM2 if they are not the order I arbitrarily like them in.
# (ADM2 is hamadryas-ness.)
if (adm[which(mds$ind_id == "97069"),]$V1 > 0.5) {
	adm = data.frame(V1=adm$V2, V2=adm$V1)
}

names(adm) = paste0("ADM", 1:2)
mds.adm = cbind(mds, adm)

# Source: http://goo.gl/K4yh
lm_eqn = function(df){
	m = lm(y ~ x, df);
	eq = substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
			list(	a = format(coef(m)[1], digits = 2), 
					b = format(coef(m)[2], digits = 2), 
					r2 = format(summary(m)$r.squared, digits = 3)))
	as.character(as.expression(eq));                 
}

# ----------------------------------------------------------------------------------------
# --- Infer "speculative" hybrid index for pure animals without one
# ----------------------------------------------------------------------------------------

mds.adm$speculative.hybrid_index = mds.adm$hybrid_index

mds.adm[mds.adm$taxon == "Hamadryas" & 
		is.na(mds.adm$hybrid_index),]$speculative.hybrid_index = 0
mds.adm[mds.adm$taxon == "Anubis"    & 
		is.na(mds.adm$hybrid_index),]$speculative.hybrid_index = 1

# ----------------------------------------------------------------------------------------
# --- Create boxplot of ADMIXTURE inferred ancestry by taxon
# ----------------------------------------------------------------------------------------

mds.adm$taxon.f = factor(mds.adm$taxon, levels = c("Anubis","Hybrid","Hamadryas"))

p = ggplot(aes(y=ADM1, x=taxon.f), data=mds.adm, pch=16) + 
	geom_boxplot() +
	xlab("Phenotypic Taxon") + ylab("ADMIXTURE Inferred Ancestry")

ggsave(p, file='results/boxplot_ancestry_by_taxon.pdf')

# ----------------------------------------------------------------------------------------
# --- Create boxplot of ADMIXTURE inferred ancestry by group
# ----------------------------------------------------------------------------------------

p = ggplot(aes(y=ADM1, x=group), data=mds.adm, pch=16) + 
	geom_boxplot() +
	xlab("Group") + ylab("ADMIXTURE Inferred Ancestry")

ggsave(p, file='results/boxplot_ancestry_by_group.pdf')

# ----------------------------------------------------------------------------------------
# --- Plot ADMIXTURE inferred ancestry as function of hybrid index
# ----------------------------------------------------------------------------------------

# Create a data frame to hold your label variables
data.label = data.frame(x=0.2, y=0.9, 
	label = c(lm_eqn(data.frame(x=mds.adm$hybrid_index, y=mds.adm$ADM1))))

p = ggplot(aes(hybrid_index, ADM1), data=mds.adm, pch=16) + 
	xlim(0,1) + ylim(0,1) + 
	xlab("Phenotypic Hybrid Index") + ylab("ADMIXTURE Inferred Ancestry") + 
	geom_point() + 
	geom_smooth(method=lm, se=FALSE) +
	geom_text(data = data.label, aes(x=x, y=y, label=label), size=4, parse=TRUE)

lm.fit = lm(mds.adm$hybrid_index ~ mds.adm$ADM1)
#summary(lm.fit)

ggsave(p, file='results/ancestries_adm_vs_phi.pdf')

# Now with inferred hybrid indices:

data.label = data.frame(x=0.2, y=0.9, 
	label = c(lm_eqn(data.frame(x=mds.adm$speculative.hybrid_index, y=mds.adm$ADM1))))

p = ggplot(aes(speculative.hybrid_index, ADM1), data=mds.adm, pch=16) + 
	xlim(0,1) + ylim(0,1) + 
	xlab("Phenotypic Hybrid Index (Inferred)") + ylab("ADMIXTURE Inferred Ancestry") + 
	geom_point() + 
	geom_smooth(method=lm, se=FALSE) +
	geom_text(data = data.label, aes(x=x, y=y, label=label), size=4, parse=TRUE)

lm.fit.spec = lm(mds.adm$speculative.hybrid_index ~ mds.adm$ADM1)
#summary(lm.fit)

ggsave(p, file='results/ancestries_adm_vs_phi_inferred.pdf')

# ----------------------------------------------------------------------------------------

# Plot residuals against phenotype-inferred ancestry and Q-Q plot
pdf(file='results/ancestries_adm_vs_phi_resid.pdf', width=8, height=4)
	layout(matrix(c(1,2),1,2))
	residualPlots(lm.fit, ~ 1, fitted=TRUE)
	qqPlot(lm.fit)
dev.off()

# Now with inferred hybrid indices:

pdf(file='results/ancestries_adm_vs_phi_resid_inferred.pdf', width=8, height=4)
	layout(matrix(c(1,2),1,2))
	residualPlots(lm.fit.spec, ~ 1, fitted=TRUE)
	qqPlot(lm.fit.spec)
dev.off()

# ----------------------------------------------------------------------------------------
# --- Plot first dimension of IBS MDS as function of hybrid index
# ----------------------------------------------------------------------------------------

# Scale C1 to 0-1 range
mds.adm$C1.norm = 1 - (mds.adm$C1-min(mds.adm$C1))/(max(mds.adm$C1)-min(mds.adm$C1))

data.label = data.frame(x=0.2, y=0.9, 
	label = c(lm_eqn(data.frame(x=mds.adm$hybrid_index, y=mds.adm$C1.norm))))

p = ggplot(aes(hybrid_index, C1.norm), data=mds.adm, pch=16) + 
	xlim(0,1) + ylim(0,1) + 
	xlab("Phenotypic Hybrid Index") + ylab("IBS MDS Dimension 1 (Scaled)") + 
	geom_point() + 
	geom_smooth(method=lm, se=FALSE) +
	geom_text(data = data.label, aes(x=x, y=y, label=label), size=4, parse=TRUE)

lm.fit = lm(mds.adm$hybrid_index ~ mds.adm$C1.norm)
#summary(lm.fit)

ggsave(p, file='results/ancestries_mds_vs_phi.pdf')

# Now with inferred hybrid indices:

data.label = data.frame(x=0.2, y=0.9, 
	label = c(lm_eqn(data.frame(x=mds.adm$speculative.hybrid_index, y=mds.adm$C1.norm))))

p = ggplot(aes(speculative.hybrid_index, C1.norm), data=mds.adm, pch=16) + 
	xlim(0,1) + ylim(0,1) + 
	xlab("Phenotypic Hybrid Index (Inferred)") + ylab("IBS MDS Dimension 1 (Scaled)") + 
	geom_point() + 
	geom_smooth(method=lm, se=FALSE) +
	geom_text(data = data.label, aes(x=x, y=y, label=label), size=4, parse=TRUE)

lm.fit.spec = lm(mds.adm$speculative.hybrid_index ~ mds.adm$C1.norm)
#summary(lm.fit)

ggsave(p, file='results/ancestries_mds_vs_phi_inferred.pdf')

# ----------------------------------------------------------------------------------------

# Plot residuals against phenotype-inferred ancestry and Q-Q plot
pdf(file='results/ancestries_mds_vs_phi_resid.pdf', width=8, height=4)
	layout(matrix(c(1,2),1,2))
	residualPlots(lm.fit, ~ 1, fitted=TRUE)
	qqPlot(lm.fit)
dev.off()

# Now with inferred hybrid indices:

pdf(file='results/ancestries_mds_vs_phi_resid_inferred.pdf', width=8, height=4)
	layout(matrix(c(1,2),1,2))
	residualPlots(lm.fit.spec, ~ 1, fitted=TRUE)
	qqPlot(lm.fit.spec)
dev.off()

# ----------------------------------------------------------------------------------------
# --- Plot first dimension of IBS MDS as function of ADMIXTURE inferred ancestry
# ----------------------------------------------------------------------------------------

data.label = data.frame(x=0.2, y=0.9, 
	label = c(lm_eqn(data.frame(x=mds.adm$ADM1, y=mds.adm$C1.norm))))

p = ggplot(aes(ADM1, C1.norm), data=mds.adm, pch=16) + 
	xlim(0,1) + ylim(0,1) + 
	xlab("ADMIXTURE Inferred Ancestry") + ylab("IBS MDS Dimension 1 (Scaled)") + 
	geom_point() + 
	geom_smooth(method=lm, se=FALSE) +
	geom_text(data = data.label, aes(x=x, y=y, label=label), size=4, parse=TRUE)

lm.fit = lm(mds.adm$ADM1 ~ mds.adm$C1.norm)
#summary(lm.fit)

ggsave(p, file='results/ancestries_mds_vs_adm.pdf')
