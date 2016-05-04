#!/bin/sh

# ------------------------------------------------------------------------------
# --- Back-up kernel smoothing script output (weighted averages files)
# ------------------------------------------------------------------------------

TODAY=`date +%Y_%m_%d`

# Compress
tar -cf - results/*weighted_averages*.txt | \
	gzip > results/pop_gen_weighted_averages_${TODAY}.tar.gz

# Back-up on AWS
module load python/intel/2.7.6
aws s3 mb s3://pop-gen-weighted-averages
aws s3 cp results/pop_gen_weighted_averages_${TODAY}.tar.gz s3://pop-gen-weighted-averages/

# Make a hard link to an archive sans date for use as dependency in Makefile
ln -f results/pop_gen_weighted_averages_${TODAY}.tar.gz results/pop_gen_weighted_averages.tar.gz
