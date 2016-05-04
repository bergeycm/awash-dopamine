#!/bin/sh

# ----------------------------------------------------------------------------------------
# --- Grab simple gene lists from find_orthologs.pl output
# ----------------------------------------------------------------------------------------

ANALYSIS=`ls results/*_weighted_averages*.orthologs.genes.out | \
          sed -e "s/results\/\(.*\)_weighted_averages.*/\1/" | sort | uniq`

# Make directory to hold gene lists
OUT_DIR=results/gene_lists/overrep_test
mkdir -p $OUT_DIR

COUNT="0"

# Combine all chromosomes, and make simple file containing, for all genes:
# weighted estimate and p-value
for i in ${ANALYSIS[*]}; do

    ls -v results/${i}_weighted_averages_chr*.orthologs.genes.out | xargs cat \
        > results/${i}.all.orthologs.genes.out
    
    echo -e "genes\tw.est\tp.val" > $OUT_DIR/${i}.all.orthologs.genes.out.est.pval
    
    # Grab just genes, weighted estimate, and p-value
    # and split lines with two or more genes (duplicating values)
    awk -v OFS='\t' '{if ($5) { print $5,$3,$4 }}' \
        results/${i}.all.orthologs.genes.out | \
        perl -lane 'foreach (split (";", $F[0])) {print "$_\t$F[1]\t$F[2]"}' \
        >> $OUT_DIR/${i}.all.orthologs.genes.out.est.pval

	# Inspiration: http://stackoverflow.com/a/16045476/1287231
	awk '{
         if ($1 == "genes") next
         sum[$1] += $2;
         if (cnt[$1]++ == 0) { max_est[$1] = $2; min_est[$1] = $2; min_p[$1] = $3 }
         if ($2 > max_est[$1]) max_est[$1] = $2
         if ($2 < min_est[$1]) min_est[$1] = $2
         if ($3 < min_p[$1]) min_p[$1] = $3
     }
     END {
         printf "%-12s %4s %10s %10s %10s %10s\n", 
         	"Gene", "N", "Min-est", "Max-est", "Avg-est", "Min-p"
         for (key in sum)
         {
             printf "%-12s %4d %1.8f %1.8f %1.8f %1.8f\n", 
             	key, cnt[key], min_est[key], max_est[key], sum[key]/cnt[key], min_p[key]
         }
     }' $OUT_DIR/${i}.all.orthologs.genes.out.est.pval \
     > $OUT_DIR/${i}.all.orthologs.genes.out.est.pval.extremes

done

exit
