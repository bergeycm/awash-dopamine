#!/bin/sh

# See https://groups.google.com/a/soe.ucsc.edu/forum/#!topic/genome/IU9QUF-dC5M

module load kent/intel
module load bedtools/intel/2.17.0

cd /scratch/cmb433/awash_pipeline-master
mkdir refGene
cd refGene

# Download refGene file
rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/rheMac3/database/refGene.txt.gz ./
gzip -d refGene.txt.gz 

# Get rid of bin column using genePredToGtf
cut -f 2- refGene.txt | genePredToGtf file stdin refGene.gtf

# Sort
sort -k1,1 -k4,4 refGene.gtf > refGene.sort.gtf
 
# Test

#	cat test.bed  
#	chr20	47473850	47579515

intersectBed -wa -a refGene.sort.gtf -b test.bed > test_genes.bed
