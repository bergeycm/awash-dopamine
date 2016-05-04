#!/bin/sh

# ------------------------------------------------------------------------------
# --- Reduce allelic depth or genotype file to just SNPs that passed filtering
# ------------------------------------------------------------------------------

IN_FILE=$1         # Has header. Allelic depth or genotype file, 
                   #   like "data/all.pass.snp.allelic_depth"
SNP_FILE=$2        # SNP map file,
                   #   like "results/strict.cleaned.map"
OUTPUT_SUFFIX=$3   # Output suffix, like "reduced" or "full"

# Grab header of AD file for use in final output file
head -n1 $IN_FILE > ${IN_FILE}_${OUTPUT_SUFFIX}

# Join first two columns (chromosome and position) of AD file 
# and sort it ignoring header
sed 's/ /_/' $IN_FILE | tail -n +2 | sort > ${IN_FILE}_sort

# Reformat first two columns of SNP file to match first columns of allelic depth file
# Note tab character (CTRL-v CTRL-i)
cut -f1,2 ${SNP_FILE} | sed -e 's/^/chr/' -e 's/	/_/' | sort > ${SNP_FILE}_sort

join ${IN_FILE}_sort ${SNP_FILE}_sort | sed 's/_/ /' >> ${IN_FILE}_${OUTPUT_SUFFIX}

rm ${IN_FILE}_sort ${SNP_FILE}_sort
