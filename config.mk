# -------------------------------------------------------------------------------------- #
# --- Configuration makefile of user-editable variables 
# -------------------------------------------------------------------------------------- #

# All paths must be absolute or relative to the awash_pipeline top directory

# -------------------------------------------------------------------------------------- #
# --- Paths to input files
# -------------------------------------------------------------------------------------- #

# The following files are expected for input:
#	${INPUT_PREFIX}.nosex
#	${INPUT_PREFIX}.map
#	${INPUT_PREFIX}.fam
#	${INPUT_PREFIX}.bim
#	${INPUT_PREFIX}.bed
#	${INPUT_PREFIX}.ped

INPUT_PREFIX=baboon.pass.snp
INPUT_PREFIX_WITHX=baboon.withX.pass.snp
INPUT_PREFIX_XONLY=baboon.chrX.pass.snp

# Path to baboon genome
GENOME_FA=genomes/papAnu2/papAnu2.fa
# Path to rhesus genome
RHESUS_GENOME_FA=genomes/rheMac3/rheMac3.fa

# Papio to rhesus liftOver chain
PAP_RHE_CHAIN=liftover_chain/papAnu2_to_rheMac3.over.chain

# Rhesus refGene GTF file
RHESUS_REFGENE=refGene/refGene.sort.gtf

# RepeatMasker output file
RM_OUT=data/papAnu2.fa.out

# Tandem Repeat Finder output file
TRF_OUT=data/papAnu2.trf.bed

# -------------------------------------------------------------------------------------- #
# --- Paths to external programs
# -------------------------------------------------------------------------------------- #

# These were updated to latest version available on HPC prior to run
# These paths are just stand-ins
R=/home/bin/r/3.0.3/intel/bin
PLINK=/home/bin/plink-1.07-x86_64
ADMIXTURE=/home/bin/admixture_linux-1.22
VCFTOOLS=/home/bin/vcftools_0.1.9/bin
TABIX=/home/bin/tabix-0.2.6
PLINKSEQ=/home/bin/plinkseq/0.08/intel/bin
BLAST=/home/bin/blast
BEDTOOLS=/home/bin/BEDTools-Version-2.13.4/bin
KENT=/home/bin/kent/x86_64
PGDSPIDER=/home/bin/PGDSpider_2.0.5.1

# -------------------------------------------------------------------------------------- #
# --- Parameters for external programs
# -------------------------------------------------------------------------------------- #

# What K to use for ADMIXTURE analysis
# (See plot of cross validation error to decide, )
ADMIX_K=2
