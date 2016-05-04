# -------------------------------------------------------------------------------------- #
# --- Makefile to run awash_pipeline. 
# -------------------------------------------------------------------------------------- #

# Get user editable variables
include config.mk

# Steps. Can be called one-by-one with something like, make download_snps
# --- filter_missing
get_missing_stats : reports/${INPUT_PREFIX}.imiss
find_repeat_regions : data/excluded_repetitive_regions.plink.txt
remove_repetitive : results/all.norep.bed
remove_missing : results/all.cleaned.bed
make_small_strict_dataset : results/strict.cleaned.ped
convert_small_dataset_to_bed : results/strict.cleaned.bed
# --- mds_plot
find_snps_in_ld : results/all.cleaned.prune.out
prune_for_ld : results/all.cleaned.LDpruned.bed
make_genome_file : results/all.cleaned.LDpruned.genome
do_mds : results/all.cleaned.LDpruned.mds
dist_matrix : results/all.cleaned.LDpruned.mdist
make_mds_plot : results/ibs_mds_plot.pdf
# --- admixture
run_admixture_cv_mode : reports/ADMIXTURE_log_cv.out
plot_admixture_cv_error : reports/ADMIXTURE_CV_plot.pdf
ADMIX_K1=$(shell expr ${ADMIX_K} \+ 1)
run_admixture_ideal_k : results/all.cleaned.LDpruned.${ADMIX_K}.Q
run_admixture_ideal_k_plus : results/all.cleaned.LDpruned.${ADMIX_K1}.Q 
run_admixture_lowest_k : results/all.cleaned.LDpruned.K_lowest_CV.Q
plot_admix_results : results/admixture_plot.pdf
# --- data_subset
get_ref_list : results/all.cleaned.just.ref
find_pure_animals : results/pure_hama_list.txt results/pure_anub_list.txt results/hybrid_list.txt
subset_pure_hama : results/pure.hamadryas.bed
subset_pure_anub : results/pure.anubis.bed
subset_hybrids : results/hybrids.bed
# --- data_conversion
convert_to_vcf : results/pure.hamadryas.vcf results/pure.anubis.vcf 
convert_to_vcf_hybrids : results/hybrids.vcf
sort_compress_pures : results/pure.hamadryas.vcf.gz.tbi results/pure.anubis.vcf.gz.tbi
sort_compress_hybrids : results/hybrids.vcf.gz.tbi
make_combined_pure_vcf : results/all.pure.vcf
# --- pop_gen_stats
calculate_pi_H : results/pure.hamadryas.sites.pi results/pure.anubis.sites.pi results/hybrids.sites.pi
explore_inbreeding : results/inbreeding_by_ancestry.pdf
calculate_Fst : results/all.pure.weir.fst
kernel_smooth : results/fst_weighted_averages_chr1.txt
pop_gen_stats_back_up : results/pop_gen_weighted_averages.tar.gz
# --- kernel_smooth_windowed
make_blast_db : ${RHESUS_GENOME_FA}.nsq
CHR_NUMS=$(shell echo {1..20})
find_orthologs :$(foreach chr,$(CHR_NUMS),results/orthologs/orthologs.chr$(chr).genes.out)
assoc_est_w_orthologs : results/fst_weighted_averages_chr1.txt.orthologs.genes.out
make_gene_list_overrep : results/gene_lists/overrep_test/fst.all.orthologs.genes.out.est.pval.extremes
make_gene_list_enrichment : results/gene_lists/enrich_test/fst.all.orthologs.genes.list.wfst
assess_ortholog_success : reports/ortholog_success.txt
est_by_gene_count : reports/fst_by_gene_count_lm.txt
# --- kernel_smooth_around_ROI_papio_genes
make_papio_refgene_bed : ${RHESUS_REFGENE}.papAnu2.bed
kernel_smooth_papio_genes : ${RHESUS_REFGENE}.papAnu2.fst.sigma50000.txt
get_ROI_gene_list : ${RHESUS_REFGENE}.papAnu2.fst.sigma150000.panther.fst
plot_smoothed_values : results/smoothed_fst_against_genome_window_vs_roi.pdf
# --- panther
panther_classify : results/panther_results/fst.all.orthologs.genes.list.panther.enrich.panther.proteinclass.txt
parse_panther_output : results/panther_results/parsed_results/windowed/overrep/panther_overrep_outliers_fst_pathway.tex
compile_panther_latex : results/panther_results/parsed_results/windowed/overrep/panther_overrep_outliers_fst_pathway.pdf
compare_methods : results/panther_results/methods_comparison/methods_comparison_fst_pathway.pdf
plot_allele_freqs : results/interesting.genes.outlier.status.csv
high_fst_introgress_outlier : results/high_fst_outlier_introgression.csv

# Group steps together
filter_missing : get_missing_stats find_repeat_regions remove_repetitive remove_missing make_small_strict_dataset convert_small_dataset_to_bed
mds_plot : find_snps_in_ld prune_for_ld make_genome_file do_mds dist_matrix make_mds_plot
admixture : run_admixture_cv_mode plot_admixture_cv_error run_admixture_ideal_k run_admixture_ideal_k_plus run_admixture_lowest_k plot_admix_results
data_subset : get_ref_list find_pure_animals subset_pure_hama subset_pure_anub subset_hybrids 
data_conversion : convert_to_vcf convert_to_vcf_hybrids sort_compress_pures sort_compress_hybrids make_combined_pure_vcf
pop_gen_stats : calculate_pi_H explore_inbreeding calculate_Fst kernel_smooth pop_gen_stats_back_up
kernel_smooth_windowed : make_blast_db find_orthologs assoc_est_w_orthologs make_gene_list_overrep make_gene_list_enrichment assess_ortholog_success est_by_gene_count
kernel_smooth_around_ROI_papio_genes : make_papio_refgene_bed kernel_smooth_papio_genes get_ROI_gene_list plot_smoothed_values
panther : panther_classify parse_panther_output compile_panther_latex compare_methods plot_allele_freqs high_fst_introgress_outlier

all : filter_missing mds_plot admixture data_subset data_conversion pop_gen_stats kernel_smooth_windowed kernel_smooth_around_ROI_papio_genes panther

SHELL_EXPORT := 

# Export Make variables to child scripts
.EXPORT_ALL_VARIABLES :

# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- Filter out missing data and variants in repetitive regions
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

# -------------------------------------------------------------------------------------- #
# --- Get stats on missingness
# -------------------------------------------------------------------------------------- #

# Stats files on missingness depend on original input file
reports/${INPUT_PREFIX}.imiss : data/${INPUT_PREFIX}.ped
	@echo "# === Getting stats on missingness ============================================ #";
	${PLINK}/plink --noweb --file data/${INPUT_PREFIX} --missing --out reports/${INPUT_PREFIX}

# -------------------------------------------------------------------------------------- #
# --- Create intervals file of repetitive regions for plink
# -------------------------------------------------------------------------------------- #

# Range file for plink depends on original input file (standing in for its map file), RM output file, and TRF output file
data/excluded_repetitive_regions.plink.txt : data/${INPUT_PREFIX}.ped ${RM_OUT} ${TRF_OUT}
	@echo "# === Creating intervals file of repetitive regions =========================== #";
	perl scripts/repeat_finder.pl 

# -------------------------------------------------------------------------------------- #
# --- Remove SNPs in repetitive regions
# -------------------------------------------------------------------------------------- #

# Binary plink without repetitive region SNPs depends on original input file and repetitive intervals file
results/all.norep.bed : data/${INPUT_PREFIX}.ped data/excluded_repetitive_regions.plink.txt
	@echo "# === Removing SNPs in repetitive regions ===================================== #";
	${PLINK}/plink --noweb --file data/${INPUT_PREFIX} --exclude data/excluded_repetitive_regions.plink.txt --range  --make-bed --out results/all.norep --missing
	cp results/all.norep.log reports/plink_repetitive_region_exclusion.log

# -------------------------------------------------------------------------------------- #
# --- Remove SNPs with too much missing data
# -------------------------------------------------------------------------------------- #

# --geno 0.1 - only include SNPs that are genotyped in at least 90% of individuals
# --mind 0.9 - only include individuals that have genotypes for at least 10% of SNPs

# Cleaned SNP binary plink (BED) file depends on original input file without repetitive regions' SNPs
results/all.cleaned.bed : results/all.norep.bed
	@echo "# === Removing SNPs with too much missing data ================================ #";
	${PLINK}/plink --noweb --bfile results/all.norep --geno 0.1 --mind 0.9 --make-bed --out results/all.cleaned --missing
	cp results/all.cleaned.log reports/plink_missing_exclusion_main.log
	# Load file to get final genotyping rate
	${PLINK}/plink --noweb --bfile results/all.cleaned --out results/all.cleaned
	cp results/all.cleaned.log reports/plink_info_full.log

# -------------------------------------------------------------------------------------- #
# --- Filter more stringently to create smaller dataset
# -------------------------------------------------------------------------------------- #

# First:
# --geno 0.05 - only include SNPs that are genotyped in at least 95% of individuals
# Then:
# --mind 0.25 - only include individuals that have genotypes for at least 75% of SNPs 
# --maf 0.1   - only include SNPs with MAF >= 0.1

# Smaller SNP plink (PED) file depends on less strictly filtered BED
results/strict.cleaned.ped : results/all.cleaned.bed
	@echo "# === Filtering stringently to create smaller dataset ========================= #";	
	${PLINK}/plink --noweb --bfile results/all.cleaned --geno 0.05 --recode --out results/strict.cleaned.tmp
	cp results/strict.cleaned.tmp.log reports/plink_missing_exclusion_strict_part1.log
	${PLINK}/plink --noweb --file results/strict.cleaned.tmp --mind 0.25 --maf 0.1 --recode --out results/strict.cleaned --missing
	cp results/strict.cleaned.log reports/plink_missing_exclusion_strict_part2.log
	rm results/strict.cleaned.tmp*

results/strict.cleaned.map : results/strict.cleaned.ped

# -------------------------------------------------------------------------------------- #
# --- Convert smaller dataset file to binary format (BED)
# -------------------------------------------------------------------------------------- #

# Smaller SNP binary plink (BED) file depends on original PED
results/strict.cleaned.bed : results/strict.cleaned.ped
	@echo "# === Converting smaller dataset SNP file to binary format ==================== #";
	${PLINK}/plink --noweb --file results/strict.cleaned --make-bed --out results/strict.cleaned
	cp results/strict.cleaned.log reports/plink_info_strict.log

# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- Multidimensional scaling plot
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

# -------------------------------------------------------------------------------------- #
# --- Find SNPs in LD to prune out
# -------------------------------------------------------------------------------------- #

# Slide a sliding window of 50 SNPs across each chromosome and calculate r^2 for each 
# pair in the window. Remove one random SNP in each pair with r^2 > 0.5.
# Parameters are window size, step, and r^2 threshold.

# List of SNPS to prune out due to LD depends on input file with missing data removed
results/all.cleaned.prune.out : results/all.cleaned.bed
	@echo "# === Finding SNPs in LD to prune out ========================================= #";
	${PLINK}/plink --noweb --bfile results/all.cleaned --indep-pairwise 50 5 0.5 --out results/all.cleaned
	cp results/all.cleaned.log reports/plink_LD_pruning_part1.log

# -------------------------------------------------------------------------------------- #
# --- Prune dataset for LD
# -------------------------------------------------------------------------------------- #

# Pruned SNP file depends on list of SNPS to prune
results/all.cleaned.LDpruned.bed : results/all.cleaned.prune.out
	@echo "# === Pruning dataset for LD ================================================== #";
	${PLINK}/plink --noweb --bfile results/all.cleaned --exclude results/all.cleaned.prune.out --make-bed --out results/all.cleaned.LDpruned
	cp results/all.cleaned.LDpruned.log reports/plink_LD_pruning_part2.log

# -------------------------------------------------------------------------------------- #
# --- Make *.genome file
# -------------------------------------------------------------------------------------- #

# Genome file depends on pruned input file
results/all.cleaned.LDpruned.genome : results/all.cleaned.LDpruned.bed
	@echo "# === Making genome file ====================================================== #";
	${PLINK}/plink --noweb --bfile results/all.cleaned.LDpruned --genome --out results/all.cleaned.LDpruned
	cp results/all.cleaned.LDpruned.log reports/plink_genome_file.log

# -------------------------------------------------------------------------------------- #
# --- Do multidimensional scaling of similarities (proportions of alleles IBS)
# -------------------------------------------------------------------------------------- #

# MDS file depends on genome file
results/all.cleaned.LDpruned.mds : results/all.cleaned.LDpruned.genome
	@echo "# === Doing multidimensional scaling of similarities ========================== #";
	${PLINK}/plink --noweb --bfile results/all.cleaned.LDpruned --read-genome results/all.cleaned.LDpruned.genome --cluster --mds-plot 2 --out results/all.cleaned.LDpruned
	cp results/all.cleaned.LDpruned.log reports/plink_mds.log

# -------------------------------------------------------------------------------------- #
# --- Generate distance matrix (1 - IBS)
# -------------------------------------------------------------------------------------- #

# Distance matrix depends on genome file
results/all.cleaned.LDpruned.mdist : results/all.cleaned.LDpruned.genome
	@echo "# === Generating distance matrix (1-IBS) ====================================== #";
	${PLINK}/plink --noweb --bfile results/all.cleaned.LDpruned --read-genome results/all.cleaned.LDpruned.genome --cluster --distance-matrix --out results/all.cleaned.LDpruned
	cp results/all.cleaned.LDpruned.log reports/plink_dist_matrix.log

# -------------------------------------------------------------------------------------- #
# --- Make MDS plot of IBS
# -------------------------------------------------------------------------------------- #

# MDS plot depends on MDS file and individual info input file
results/ibs_mds_plot.pdf : results/all.cleaned.LDpruned.mds data/individual_info.csv
	@echo "# === Creating MDS plot of IBS ================================================ #";
	${R}/Rscript scripts/plot_ibs_mds.R

# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- ADMIXTURE
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

# -------------------------------------------------------------------------------------- #
# --- Run ADMIXTURE in cross-validation mode
# -------------------------------------------------------------------------------------- #

# ADMIXTURE cross-validation log depends on pruned input file
reports/ADMIXTURE_log_cv.out : results/all.cleaned.LDpruned.bed
	@echo "# === Running ADMIXTURE in cross-validation mode ============================== #";
	sh scripts/run_admixture_cv.sh

# -------------------------------------------------------------------------------------- #
# --- Plot cross-validation error as function of K (find ideal K via visual inspection)
# -------------------------------------------------------------------------------------- #

# ADMIXTURE CV plot depends on ADMIXTURE cross-validation log
reports/ADMIXTURE_CV_plot.pdf : reports/ADMIXTURE_log_cv.out
	@echo "# === Creating plot of ADMIXTURE cross-validation error ======================= #";
	${R}/Rscript scripts/plot_admixture_cv_error.R

# -------------------------------------------------------------------------------------- #
# --- Run ADMIXTURE with chosen value of K and K+1 and K with lowest CV
# -------------------------------------------------------------------------------------- #

# Add -j[num_threads] to run in parallel.

# ADMIXTURE output depends on pruned input file
results/all.cleaned.LDpruned.${ADMIX_K}.Q : results/all.cleaned.LDpruned.bed
	@echo "# === Running ADMIXTURE for chosen K ========================================== #";
	${ADMIXTURE}/admixture -B200 results/all.cleaned.LDpruned.bed ${ADMIX_K}
	mv all.cleaned.LDpruned.${ADMIX_K}.Q results/
	mv all.cleaned.LDpruned.${ADMIX_K}.Q_se results/
	mv all.cleaned.LDpruned.${ADMIX_K}.Q_bias results/
	mv all.cleaned.LDpruned.${ADMIX_K}.P results/

results/all.cleaned.LDpruned.${ADMIX_K1}.Q : results/all.cleaned.LDpruned.bed
	@echo "# === Running ADMIXTURE for chosen K, plus 1 ================================== #";
	${ADMIXTURE}/admixture -B200 results/all.cleaned.LDpruned.bed ${ADMIX_K1}
	mv all.cleaned.LDpruned.${ADMIX_K1}.Q results/
	mv all.cleaned.LDpruned.${ADMIX_K1}.Q_se results/
	mv all.cleaned.LDpruned.${ADMIX_K1}.Q_bias results/
	mv all.cleaned.LDpruned.${ADMIX_K1}.P results/

ADMIX_LOW_CV_K:=$(shell awk 'BEGIN { FS = ":" } ; NR == 1 || $$2 < min {line = $$0; min = $$2}END{print line}' reports/ADMIXTURE_log_cv.out | sed -e "s/.*K=\([0-9]*\).*/\1/")
results/all.cleaned.LDpruned.K_lowest_CV.Q : results/all.cleaned.LDpruned.bed
	@echo "# === Running ADMIXTURE for K with lowest CV ================================== #";
	${ADMIXTURE}/admixture -B200 results/all.cleaned.LDpruned.bed ${ADMIX_LOW_CV_K}
	mv all.cleaned.LDpruned.${ADMIX_LOW_CV_K}.Q      results/
	mv all.cleaned.LDpruned.${ADMIX_LOW_CV_K}.Q_se   results/
	mv all.cleaned.LDpruned.${ADMIX_LOW_CV_K}.Q_bias results/
	mv all.cleaned.LDpruned.${ADMIX_LOW_CV_K}.P      results/
	cd results
	ln -s all.cleaned.LDpruned.${ADMIX_LOW_CV_K}.Q      ./all.cleaned.LDpruned.K_lowest_CV.Q
	ln -s all.cleaned.LDpruned.${ADMIX_LOW_CV_K}.Q_se   ./all.cleaned.LDpruned.K_lowest_CV.Q_se
	ln -s all.cleaned.LDpruned.${ADMIX_LOW_CV_K}.Q_bias ./all.cleaned.LDpruned.K_lowest_CV.Q_bias
	ln -s all.cleaned.LDpruned.${ADMIX_LOW_CV_K}.P      ./all.cleaned.LDpruned.K_lowest_CV.P
	cd ..

# -------------------------------------------------------------------------------------- #
# --- Plot ADMIXTURE results
# -------------------------------------------------------------------------------------- #

# ADMIXTURE plot depends on ADMIXTURE output file
results/admixture_plot.pdf : results/all.cleaned.LDpruned.${ADMIX_K}.Q
	@echo "# === Plotting ADMIXTURE results ============================================== #";
	${R}/Rscript scripts/plot_admixture_results.R

# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- Data slicing and conversion
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

# -------------------------------------------------------------------------------------- #
# --- Figure out which allele plink considers the "reference"
# -------------------------------------------------------------------------------------- #

# List of reference alleles depends on pruned binary PED file
results/all.cleaned.just.ref : results/all.cleaned.bed
	@echo "# === Making list of reference alleles ======================================== #";
	${PLINK}/plink --noweb --bfile results/all.cleaned --recode-lgen --with-reference --out results/all.cleaned
	cp results/all.cleaned.log reports/plink_ref_alleles.log
	cut -d' ' -f1,2 results/all.cleaned.ref > results/all.cleaned.just.ref

# -------------------------------------------------------------------------------------- #
# --- Find pure anubis and hamadryas (and get list of hybrids too)
# -------------------------------------------------------------------------------------- #

# Lists of pure animals depend on ADMIXTURE output file
results/%_list.txt : results/all.cleaned.LDpruned.${ADMIX_K}.Q
	@echo "# === Finding pure animals ==================================================== #";
	${R}/Rscript scripts/find_pure_animals.R

results/pure_anub_list.txt : results/pure_hama_list.txt
results/hybrid_list.txt : results/pure_hama_list.txt

# -------------------------------------------------------------------------------------- #
# --- Make subsets of just pure anubis and hamadryas (forcing ref allele to be same)
# -------------------------------------------------------------------------------------- #

# Pure Hamadryas subset depends on list of pure Hamadryas animals
results/pure.hamadryas.bed : results/pure_hama_list.txt
	@echo "# === Pulling out pure hamadryas ============================================== #";
	${PLINK}/plink --noweb --bfile results/all.cleaned --keep results/pure_hama_list.txt --reference-allele results/all.cleaned.just.ref --make-bed --out results/pure.hamadryas
	cp results/pure.hamadryas.log reports/plink_subset_hama.log

# Pure Anubis subset depends on list of pure Anubis animals
results/pure.anubis.bed : results/pure_anub_list.txt
	@echo "# === Pulling out pure anubis ================================================= #";
	${PLINK}/plink --noweb --bfile results/all.cleaned --keep results/pure_anub_list.txt --reference-allele results/all.cleaned.just.ref --make-bed --out results/pure.anubis
	cp results/pure.anubis.log reports/plink_subset_anub.log

# -------------------------------------------------------------------------------------- #
# --- Make subsets of just hybrids (forcing ref allele to be same)
# -------------------------------------------------------------------------------------- #

# Hybrids subset depends on list of hybrid animals
results/hybrids.bed : results/hybrid_list.txt
	@echo "# === Pulling out hybrids ===================================================== #";
	${PLINK}/plink --noweb --bfile results/all.cleaned --keep results/hybrid_list.txt --reference-allele results/all.cleaned.just.ref --make-bed --out results/hybrids
	cp results/hybrids.log reports/plink_subset_hybr.log

# -------------------------------------------------------------------------------------- #
# --- Convert to VCF format (pure anubis, pure hamadryas)
# -------------------------------------------------------------------------------------- #

# Pure VCF file depends on pure BED file
results/pure.%.vcf : results/pure.%.bed
	@echo "# === Converting pure animals' PED file to VCF format ========================= #";
	${PLINKSEQ}/pseq $(basename $@).pseq new-project
	${PLINKSEQ}/pseq $(basename $@).pseq load-plink --file $(basename $@) --id $(basename $@)
	${PLINKSEQ}/pseq $(basename $@).pseq write-vcf > $@
	rm -r $(basename $@).pseq*

# -------------------------------------------------------------------------------------- #
# --- Convert to VCF format (hybrids)
# -------------------------------------------------------------------------------------- #

# Pure VCF file depends on pure BED file
results/hybrids.vcf : results/hybrids.bed
	@echo "# === Converting hybrid animals' PED file to VCF format ======================= #";
	${PLINKSEQ}/pseq $(basename $@).pseq new-project
	${PLINKSEQ}/pseq $(basename $@).pseq load-plink --file $(basename $@) --id $(basename $@)
	${PLINKSEQ}/pseq $(basename $@).pseq write-vcf > $@
	rm -r $(basename $@).pseq*

# -------------------------------------------------------------------------------------- #
# --- Sort and compress pure VCF files
# -------------------------------------------------------------------------------------- #

# Sorted and compressed VCF files depend on original VCF files
results/pure.%.vcf.gz.tbi : results/pure.%.vcf
	@echo "# === Sorting and compressing pure VCF files ================================== #";
	${VCFTOOLS}/vcf-sort $< > $<.sort
	${TABIX}/bgzip -c $<.sort > $(basename $@)
	rm $<.sort
	${TABIX}/tabix -p vcf $(basename $@)

# -------------------------------------------------------------------------------------- #
# --- Sort and compress hybrid VCF files
# -------------------------------------------------------------------------------------- #

# Sorted and compressed VCF files depend on original VCF files
results/hybrids.vcf.gz.tbi : results/hybrids.vcf
	@echo "# === Sorting and compressing hybrid VCF files ================================ #";
	${VCFTOOLS}/vcf-sort $< > $<.sort
	${TABIX}/bgzip -c $<.sort > $(basename $@)
	rm $<.sort
	${TABIX}/tabix -p vcf $(basename $@)

# -------------------------------------------------------------------------------------- #
# --- Make file containing all pure animals
# -------------------------------------------------------------------------------------- #

PATH := ${PATH}:${TABIX}

# VCF file with all pure animals depends on separate compressed VCF files of pure animals
results/all.pure.vcf : results/pure.hamadryas.vcf.gz.tbi results/pure.anubis.vcf.gz.tbi
	@echo "# === Combining pure VCF files ================================================ #";
	${VCFTOOLS}/vcf-merge $(basename $^) > $@

# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- Pop gen stats calculation
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

# -------------------------------------------------------------------------------------- #
# --- Calculate pi, nucleotide diversity, and H, observed heterozygosity, among others
# -------------------------------------------------------------------------------------- #

#	--freq              per-site frequency              *.frq
#	--counts            per-site counts                 *.frq.count
#	--depth             mean depth per individual       *.idepth            (fails)
#	--site-depth        depth for each site (summed)    *.ldepth            (fails)
#	--site-mean-depth   depth for each site (mean)      *.ldepth.mean       (fails)
#	--geno-depth        depth for each genotype         *.gdepth            (fails)
#	--site-quality      per-site SNP quality            *.lqual             (fails)
#	--het               inbreeding coefficient, F       *.het
#	--hardy             p-value for HWE                 *.hwe
#	--missing           missingness per ind, per site   *.lmiss, out.imiss
#	--geno-r2           squard corr coef btwn genotypes *.geno.ld           (takes awhile)
#	--SNPdensity <int>  SNP # and density for bin size  *.snpden
#	--TsTv <int>        Ts / Tv ratio for bin size      *.TsTv, out.TsTv.summary
#	--singletons        location of singletons          *.singletons
#	--site-pi           nuc. diversity per-site         *.sites.pi
#	--window-pi <int>   nuc. diversity in windows       *.windowed.pi
#	--TajimaD <int>     TajD in bins of size <int>      *.Tajima.D

# Output pop gen stats files (using *.pi as stand-in for all) depend on fixed VCF file
# Pures:
results/pure.%.sites.pi : results/pure.%.vcf
	@echo "# === Calculating nucleotide diversity and heterozygosity (pures) ============= #";
	${VCFTOOLS}/vcftools --vcf $< --out $(basename $(basename $@)) --freq
	${VCFTOOLS}/vcftools --vcf $< --out $(basename $(basename $@)) --counts
	${VCFTOOLS}/vcftools --vcf $< --out $(basename $(basename $@)) --het
	${VCFTOOLS}/vcftools --vcf $< --out $(basename $(basename $@)) --hardy
	${VCFTOOLS}/vcftools --vcf $< --out $(basename $(basename $@)) --missing
	${VCFTOOLS}/vcftools --vcf $< --out $(basename $(basename $@)) --SNPdensity 1000000
	${VCFTOOLS}/vcftools --vcf $< --out $(basename $(basename $@)) --TsTv 1000000
	${VCFTOOLS}/vcftools --vcf $< --out $(basename $(basename $@)) --singletons
	${VCFTOOLS}/vcftools --vcf $< --out $(basename $(basename $@)) --site-pi
	${VCFTOOLS}/vcftools --vcf $< --out $(basename $(basename $@)) --TajimaD 300000

results/pure.anubis.frq.count : results/pure.anubis.sites.pi
results/pure.hamadryas.frq.count : results/pure.hamadryas.sites.pi

# Hybrids:
results/hybrids.sites.pi : results/hybrids.vcf
	@echo "# === Calculating nucleotide diversity and heterozygosity (hybrids) =========== #";
	${VCFTOOLS}/vcftools --vcf $< --out $(basename $(basename $@)) --freq
	${VCFTOOLS}/vcftools --vcf $< --out $(basename $(basename $@)) --counts
	${VCFTOOLS}/vcftools --vcf $< --out $(basename $(basename $@)) --het
	${VCFTOOLS}/vcftools --vcf $< --out $(basename $(basename $@)) --hardy
	${VCFTOOLS}/vcftools --vcf $< --out $(basename $(basename $@)) --missing
	${VCFTOOLS}/vcftools --vcf $< --out $(basename $(basename $@)) --SNPdensity 1000000
	${VCFTOOLS}/vcftools --vcf $< --out $(basename $(basename $@)) --TsTv 1000000
	${VCFTOOLS}/vcftools --vcf $< --out $(basename $(basename $@)) --singletons
	${VCFTOOLS}/vcftools --vcf $< --out $(basename $(basename $@)) --site-pi
	${VCFTOOLS}/vcftools --vcf $< --out $(basename $(basename $@)) --TajimaD 300000

results/hybrids.frq.count : results/hybrids.sites.pi

# ----------------------------------------------------------------------------------------
# --- Explore inbreeding levels
# ----------------------------------------------------------------------------------------

# Last inbreeding plot depends on hybrids' diversity output file (stand-in for *.het files)
results/inbreeding_by_ancestry.pdf : results/hybrids.sites.pi
	@echo "# === Exploring inbreeding levels ============================================= #";
	${R}/Rscript scripts/explore_inbreeding.R

# -------------------------------------------------------------------------------------- #
# --- Calculate Fst
# -------------------------------------------------------------------------------------- #

# Need lists of animals in each population
# This generates out.pures.weir.fst
# Based on Weir and Cockerham 1984

# Output Fst file depends on population lists and combined pures VCF file
results/all.pure.weir.fst : results/all.pure.vcf results/pure_hama_list.txt results/pure_anub_list.txt
	@echo "# === Calculating Fst based on Weir and Cockerham 1984 ======================== #";
	${VCFTOOLS}/vcftools --vcf results/all.pure.vcf --out results/all.pure --weir-fst-pop results/pure_hama_list.txt --weir-fst-pop results/pure_anub_list.txt

# -------------------------------------------------------------------------------------- #
# --- Do kernel smoothing of population genetic stats (currently skipped)
# -------------------------------------------------------------------------------------- #

CHRS=$(shell echo chr{1..20})

# Output Fst weighted averages file (standing in for pi and het) depends on results of windowed kernel smoothing
results/fst_weighted_averages_chr1.txt : results/pure.hamadryas.sites.pi results/pure.anubis.sites.pi results/pure.hamadryas.hwe results/pure.anubis.hwe results/all.pure.weir.fst
	@echo "# === Smoothing the population genetic stats ================================== #";
	${R}/Rscript scripts/kernel_smooth.R ${CHRS} --pi_h --pi_a --het_h --het_a --fst

# -------------------------------------------------------------------------------------- #
# --- Back-up kernel smoothing script output (weighted averages files) (currently skipped)
# -------------------------------------------------------------------------------------- #

# Compressed archive of weighted averages results files depends on one output weighted averages file
results/pop_gen_weighted_averages.tar.gz : results/fst_weighted_averages_chr1.txt
	@echo "# === Backing up weighted average results ===================================== #";
	sh scripts/backup_weighted_avgs.sh

# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- Find genes near interesting SNPs
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

# -------------------------------------------------------------------------------------- #
# --- Make BLAST database for rhesus genome
# -------------------------------------------------------------------------------------- #

# BLASTable database depends on rhesus genome FASTA
${RHESUS_GENOME_FA}.nsq : ${RHESUS_GENOME_FA}
	@echo "# === Making BLAST database for rhesus genome ================================= #";
	${BLAST}/formatdb -p F -i ${RHESUS_GENOME_FA} -n ${RHESUS_GENOME_FA}

# -------------------------------------------------------------------------------------- #
# --- Find orthologs of windows (Only needs to be called on one estimate, say Fst.)
# -------------------------------------------------------------------------------------- #

# List of orthologous regions and genes contained therein depends on the Fst measures (stand-in), refGene file, chain file, and scripts/find_orthologs.pl
results/orthologs/orthologs.chr%.genes.out : results/fst_weighted_averages_chr%.txt ${PAP_RHE_CHAIN} ${RHESUS_REFGENE} #scripts/find_orthologs.pl
	@echo "# === Finding orthologs of windows - chr $* =================================== #";
	perl scripts/find_orthologs.pl $<

# -------------------------------------------------------------------------------------- #
# --- Associate smoothed estimates with their nearby genes
# -------------------------------------------------------------------------------------- #

# File of smoothed estimates, p-values, and nearby genes (chr1 as stand-in) depends on the smoothed Fst (stand-in) and file of orthologs
results/fst_weighted_averages_chr1.txt.orthologs.genes.out : results/orthologs/orthologs.chr1.genes.out results/fst_weighted_averages_chr1.txt
	@echo "# === Associating smoothed estimates with their nearby genes ================== #";
	perl scripts/assoc_ests_w_orthologs.pl

# -------------------------------------------------------------------------------------- #
# --- Get simple lists of genes, weighted estimates, and p-values and of extreme values
# -------------------------------------------------------------------------------------- #

# Gene list (one of two simple gene lists output) depends on *.genes.out files from find_orthologs.pl
results/gene_lists/overrep_test/fst.all.orthologs.genes.out.est.pval.extremes : results/fst_weighted_averages_chr1.txt.orthologs.genes.out
	@echo "# === Making simple gene list ================================================= #";
	sh scripts/get_gene_list.sh

results/fst.all.orthologs.genes.out : results/gene_lists/overrep_test/fst.all.orthologs.genes.out.est.pval.extremes

# After this step, run statistical overrepresentation tests manually in PANTHER.
# See scripts/steps_panther_fst_outliers.md

# -------------------------------------------------------------------------------------- #
# --- Get gene list with smoothed Fst for enrichment analysis
# -------------------------------------------------------------------------------------- #

# Gene list with Fst depends on results/fst.all.orthologs.genes.out (output during get_gene_list.sh)
results/gene_lists/enrich_test/fst.all.orthologs.genes.list.wfst : results/fst.all.orthologs.genes.out
	@echo "# === Outputting gene list with average of smoothed values ==================== #";
	perl scripts/get_gene_list_for_enrichment.pl

# -------------------------------------------------------------------------------------- #
# --- Assess how well ortholog finding worked
# -------------------------------------------------------------------------------------- #

# Report depends on output file from get_gene_list.sh (stands in for other files)
reports/ortholog_success.txt : results/gene_lists/overrep_test/fst.all.orthologs.genes.out.est.pval.extremes
	@echo "# === Assessing ortholog finding success ====================================== #";
	sh scripts/assess_ortholog_success.sh

# -------------------------------------------------------------------------------------- #
# --- Look at how smoothed Fst changes by region's gene count
# -------------------------------------------------------------------------------------- #

# Report depends on results/fst.all.orthologs.genes.out (output during get_gene_list.sh)
reports/fst_by_gene_count_lm.txt : results/fst.all.orthologs.genes.out
	@echo "# === Exploring how estimates change with gene count ========================== #";
	${R}/Rscript scripts/est_by_gene_count.R

# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- Do kernel smoothing of Fst around refGene genes
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

# -------------------------------------------------------------------------------------- #
# --- Convert Macaque refGene GTF file to baboon BED file
# -------------------------------------------------------------------------------------- #

# Papio refGene BED file depends on 
${RHESUS_REFGENE}.papAnu2.bed : ${RHESUS_REFGENE} liftover_chain/rheMac3_to_papAnu2.over.chain #scripts/get_bed_of_genes_from_gtf.pl 
	@echo "# === Convert Macaque refGene GTF file to baboon BED file ===================== #";
	perl scripts/get_bed_of_genes_from_gtf.pl ${RHESUS_REFGENE} > ${RHESUS_REFGENE}.bed
	${KENT}/liftOver ${RHESUS_REFGENE}.bed \
		liftover_chain/rheMac3_to_papAnu2.over.chain \
		${RHESUS_REFGENE}.papAnu2.bed.tmp \
		${RHESUS_REFGENE}.papAnu2.unmapped.bed \
		-minMatch=0.5
	echo -e "Chr\tBegin\tEnd\tGeneID" > ${RHESUS_REFGENE}.papAnu2.bed
	cat ${RHESUS_REFGENE}.papAnu2.bed.tmp >> ${RHESUS_REFGENE}.papAnu2.bed

# -------------------------------------------------------------------------------------- #
# --- Kernel smooth around ROI: Papio genes
# -------------------------------------------------------------------------------------- #

# Results file of kernel smoothed Fst estimates depends on baboon refGene BED file and original Fst estimates
${RHESUS_REFGENE}.papAnu2.fst.sigma50000.txt : ${RHESUS_REFGENE}.papAnu2.bed results/all.pure.weir.fst #scripts/kernel_smooth_around_ROI.R
	@echo "# === Kernel smooth around ROI: Papio genes =================================== #";
	echo "Kernel smoothing around ROI must be run in parallel outside of the Makefile." | tee /dev/stderr
	cat scripts/cmds_to_submit_KS_ROI_jobs.txt | tee /dev/stderr
	exit 2

# -------------------------------------------------------------------------------------- #
# --- Make PANTHER enrichment test input file from estimates smoothed around genes
# -------------------------------------------------------------------------------------- #

# Results file of kernel smoothed estimates (Fst stand-in) depends on baboon refGene BED file and original estimates
${RHESUS_REFGENE}.papAnu2.fst.sigma150000.panther.fst : ${RHESUS_REFGENE}.papAnu2.fst.sigma50000.txt
	@echo "# === Making PANTHER enrichment test input from ests. smoothed around genes === #";
	${R}/Rscript scripts/get_gene_list_ROI.R "refGene/refGene.sort.gtf.papAnu2.bed" "fst"

# -------------------------------------------------------------------------------------- #
# --- Plot smoothed values
# -------------------------------------------------------------------------------------- #

# Plot of smoothed estimates across genome (Fst as stand-in) depends on results file of kernel smoothed estimates
results/smoothed_fst_against_genome_window_vs_roi.pdf : ${RHESUS_REFGENE}.papAnu2.fst.sigma150000.panther.fst
	@echo "# === Plotting smoothed values ================================================ #";
	${R}/Rscript scripts/plot_smoothed_values.R

# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- Classify genes with PANTHER and do tests for interesting annotations
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

# -------------------------------------------------------------------------------------- #
# --- Classify genes with PANTHER clone and do enrichment test and overrep test
# -------------------------------------------------------------------------------------- #

# PANTHER output file (place holder) depends on gene list w/ smoothed values for window and ROI
results/panther_results/fst.all.orthologs.genes.list.panther.enrich.panther.proteinclass.txt : results/gene_lists/enrich_test/fst.all.orthologs.genes.list.wfst ${RHESUS_REFGENE}.papAnu2.fst.sigma150000.panther.fst
	@echo "# === Classifying genes with PANTHER and doing enrichment test ================ #";
	${R}/Rscript scripts/panther_classify.R

# -------------------------------------------------------------------------------------- #
# --- Parse PANTHER enrichment and overrepresentation test output
# -------------------------------------------------------------------------------------- #

# LaTeX file created from PANTHER test results depends on PANTHER output
results/panther_results/parsed_results/windowed/overrep/panther_overrep_outliers_fst_pathway.tex : results/panther_results/fst.all.orthologs.genes.list.panther.enrich.panther.proteinclass.txt
	@echo "# === Parsing PANTHER test output ============================================= #";
	${R}/Rscript scripts/parse_panther_results.R

# -------------------------------------------------------------------------------------- #
# --- Compile PDFs summarizing PANTHER results
# -------------------------------------------------------------------------------------- #

# PDF depends on LaTeX file created from PANTHER test results
results/panther_results/parsed_results/windowed/overrep/panther_overrep_outliers_fst_pathway.pdf : results/panther_results/parsed_results/windowed/overrep/panther_overrep_outliers_fst_pathway.tex
	@echo "# === Compiling PDFs summarizing PANTHER results ============================== #";
	sh scripts/build_panther_tables.sh

# -------------------------------------------------------------------------------------- #
# --- Compare PANTHER results from different smoothing methods
# -------------------------------------------------------------------------------------- #

# Output PDF (stand-in) depends on output of PANTHER (stand-in)
results/panther_results/methods_comparison/methods_comparison_fst_pathway.pdf : results/panther_results/fst.all.orthologs.genes.list.panther.enrich.panther.proteinclass.txt
	@echo "# === Comparing PANTHER results from different smoothing methods ============== #";
	Rscript scripts/compare_panther_results.R

# -------------------------------------------------------------------------------------- #
# --- Plot allele frequencies of interesting genes
# -------------------------------------------------------------------------------------- #

# Output text file depends on output of PANTHER (stand-in) and list of interesting genes
results/interesting.genes.outlier.status.csv : results/panther_results/fst.all.orthologs.genes.list.panther.enrich.panther.proteinclass.txt data/interesting_pathways_genes.txt
	@echo "# === Plotting allele frequencies of interesting genes ======================== #";
	perl scripts/plot_allele_freqs_interesting_genes.pl
	Rscript scripts/output.very.interesting.genes.R

# -------------------------------------------------------------------------------------- #
# --- Find differentiated (high Fst) introgression outliers
# -------------------------------------------------------------------------------------- #

# Output file depends on results file of kernel smoothed estimates (stand-in)
results/high_fst_outlier_introgression.csv : ${RHESUS_REFGENE}.papAnu2.fst.sigma150000.panther.fst
	@echo "# === Finding differentiated (high Fst) introgression outliers ================ #";
	Rscript scripts/find_differentiated_introgression_outliers.R
