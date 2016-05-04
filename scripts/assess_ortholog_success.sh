#!/bin/sh

# ========================================================================================
# --- For how many regions were orthologs successfully found? 
# ========================================================================================

# Print header
OUT="reports/ortholog_success.txt"
echo -e "SET\tSUCCESSES\tATTEMPTS\tPERCENT" > $OUT

# Total successes:
TOT_SUCCESS=`wc -l results/fst.all.orthologs.genes.out| cut -d' ' -f1`

# Total attempted:
TOT_ATTEMPT=`cat results/fst_weighted_averages_chr*.txt | wc -l`

# Percentage successful:
TOT_PERC=`printf '%f\n' $(echo "scale = 10; $TOT_SUCCESS / $TOT_ATTEMPT" | bc)`

# Print results for all chromosomes

echo -e "Total\t$TOT_SUCCESS\t$TOT_ATTEMPT\t$TOT_PERC" >> $OUT

# Now split by chromosome. 

for f in $( ls -v results/fst_weighted_averages_chr*.txt ); do
	
	# Figure out chromosome ID
	CHR=`echo $f | sed "s/.*\(chr[0-9]\+\).*/\1/"`
	
	# Successes for chromosome:
	CHR_SUCCESS=`wc -l ${f}.orthologs.genes.out | cut -d' ' -f1`

	# Attempts for chromosome:
	CHR_ATTEMPT=$((`wc -l $f | cut -d' ' -f1` - 1))

	# Percentage successful:
	CHR_PERC=`printf '%f\n' $(echo "scale = 10; $CHR_SUCCESS / $CHR_ATTEMPT" | bc)`

	echo -e "$CHR\t$CHR_SUCCESS\t$CHR_ATTEMPT\t$CHR_PERC" >> $OUT

done
        