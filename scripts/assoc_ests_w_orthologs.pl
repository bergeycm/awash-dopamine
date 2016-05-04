#!/usr/bin/perl

use strict;
use warnings;

use List::MoreUtils qw(uniq);

use Data::Dumper;

# ========================================================================================
# --- Script to combine genes found in macaque orthologs and smoothed estimates
# ========================================================================================

# Get all smoothed estimates files

my @smoothed_est_files = glob("results/*_weighted_averages_chr*.txt");

# Figure out the analyses and chromosomes represented from the filenames

my @analyses = 	map {s:results/(.*)_weighted_averages_chr.*\.txt:$1:g; $_;} 
				map {$_} @smoothed_est_files;
@analyses = uniq @analyses;

my @chroms = 	map {s:results/.*_weighted_averages_chr(.*)\.txt:$1:g; $_;} 
				map {$_} @smoothed_est_files;
@chroms = sort {$a <=> $b} (uniq @chroms);

# Combine smoothed estimates and orthologs

foreach my $analysis (@analyses) {

	foreach my $chr (@chroms) {
	
		print STDERR "--- Processing estimates of [$analysis] for chr$chr ---\n";
	
		my $est_file   = "results/${analysis}_weighted_averages_chr${chr}.txt";
		my $ortho_file = "results/orthologs/orthologs.chr${chr}.genes.out";
		
		# Figure out and open combined output file
		my $out_file = $est_file . ".orthologs.genes.out";
		
		open (OUT, ">$out_file")
			or die "ERROR: Could not open output file, $out_file. $!\n";
		
		# Read in file of smoothed estimates
		open (ESTS, "<$est_file")
			or die "ERROR: Could not open file of smoothed estimates, $est_file. $!\n";
		
		# Put estimate info into a hash. If line is:
		#  chr20 550001 0.160324000673071 0.8199
		# Then:
		#  %estimates{"550001"} = "0.160324000673071\t0.8199"
		
		my %estimates;
		my $header = <ESTS>;
		foreach (<ESTS>) {
		
			if (/chr\d+\s(\d+)\s([\d\.\-NA]+\s[\d\.NA]+)/) {
				$estimates{$1} = $2;
			} else {
				warn "ERROR: Malformed line in smoothed estimates! [$_]\n";
			}
		}
		
		# Now read through ortholog file, adding in estimate info 
		# and outputting to a new file.
		# Read in file of smoothed estimates
		open (ORTHO, "<$ortho_file")
			or die "ERROR: Could not open file of orthologs, $ortho_file. $!\n";
		
		# Line example:
		#  chr1:71500001-71600001	chr1:72828900-72926103	
		# Or:
		#  chr1:67600001-67700001	chr1:68900161-69000738	LEPR;

		foreach (<ORTHO>) {
		
			chomp;
		
			if (/	
					(chr\d+:(\d+)-(\d+))
					\t
					(chr\d+:\d+-\d+)
					\t
					(.*)
				/x) {
				my $baboon       = $1;
				my $baboon_start = $2;
				my $baboon_end   = $3;
				my $rhesus       = $4;
				my $gene_str     = $5;
				
				my $midpoint = ($baboon_end + $baboon_start) / 2;
				
				print STDERR "Ortholog: [$_]\n";
				print STDERR "\tBaboon: $baboon\n";
				print STDERR "\t- Baboon Start: $baboon_start\n";
				print STDERR "\t- Baboon End:   $baboon_end\n";
				print STDERR "\t- Midpoint:     $midpoint\n";
				print STDERR "\tRhesus: $rhesus\n";
				print STDERR "\tGenes: [$gene_str]\n\n";
				
				# Find match in smoothed estimates hash (weighted avg and p-val)
				my $est_info = $estimates{"$midpoint"};
				$est_info =~ s/ /\t/g;
				print STDERR "\tEstimate avg and p-val:\t[$est_info]\n\n";
				
				# Output file
				print OUT "$baboon\t$rhesus\t$est_info\t$gene_str\n";
				
			} else {
				warn "ERROR: Malformed line in orthologs file! [$_]\n";
			}
		}
		
		close OUT;
	}
}

exit;
