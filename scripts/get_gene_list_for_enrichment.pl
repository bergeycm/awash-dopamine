#!/usr/bin/perl

use strict;
use warnings;

use List::Util qw(sum);
use List::MoreUtils qw(uniq);
use Getopt::Long;

my $verbose = '';
GetOptions ('verbose' => \$verbose);

# Make output directory
my $outdir = "results/gene_lists/enrich_test";
unless(-e $outdir or mkdir $outdir) {
	die "Unable to create output directory, [$outdir]\n";
}

# Figure out analyses for which to make enrichment gene lists
# Using the previously-created simple gene lists
my @gene_lists = glob("results/*.all.orthologs.genes.out");

# Figure out the analyses represented from the filenames

my @analyses = 	map {s:results/(.*).all.orthologs.genes.out:$1:g; $_;} 
				map {$_} @gene_lists;
@analyses = uniq @analyses;

for my $analysis (@analyses) {
	
	my $est_smoothed = "results/${analysis}.all.orthologs.genes.out";
	
	#	chr1:3200001-3300001	chr1:3121054-3214747	0.633369055993775	0.000000899999910000009	CPLX1;
	#	chr1:3300001-3400001	chr1:3214847-3318280	0.539590708639337	0.0000175999982400002	CPLX1;PCGF3;
	#	chr1:3400001-3500001	chr1:3321054-3423460	0.467782341235794	0.000242899975710002	
	#	chr1:3600001-3700001	chr1:3561132-3654705	0.380691828470087	0.00723359927664007
	
	# Loop through file, putting data in temporary hash:
	# - Key is gene code, value is array of Fst values
	# Ultimately hash should be in the following format:
	# - Key is gene code, value is (average) Fst p-value.
	# Genes that fall in multiple regions have to be corrected to average p-value
	
	my %gene_est_tmp;
	my %gene_est;
	
	open (EST, "<$est_smoothed")
		or die "ERROR: Could not open file of smoothed estimates and p-values. $!\n";
	
	while (<EST>) {
	
		if (/chr[^\s]+\tchr[^\s]+\t([\d\.]+)\t([\d\.]+)\t(.*)$/) {

			my $est           = $1;
			my $pval          = $2;
			my $gene_code_str = $3;
			
			if ($gene_code_str) {
			
				print STDERR "$gene_code_str has est $est with p=$pval\n" if $verbose;
				
				# Split gene codes
				my @gene_codes = split(';', $gene_code_str);
				
				# If gene code is not yet in hash, add estimate as one-element array
				# otherwise, tack this estimate onto end of array.
				foreach (@gene_codes) {
					push(@{$gene_est_tmp{$_}}, $est);
					
				}
			}
		}
	}
	
	close EST;
	
	# Now loop through temp. hash and compute avg. estimate for genes in multiple regions
	
	foreach my $gene (keys %gene_est_tmp) {
		
		my @these_ests = @{$gene_est_tmp{$gene}};
		
		my $avg_est = sum(@these_ests) / @these_ests;
		print STDERR "\tGene $gene has avg. estimate of:\t" . $avg_est . "\n" if $verbose;
		
		# Add result to final hash
		$gene_est{$gene} = $avg_est;
	
	}
	
	# Loop through final hash and output results to file
	my $main_out = "${outdir}/${analysis}.all.orthologs.genes.list.w${analysis}";
	
	open (MAINOUT, ">$main_out")
		or die "ERROR: Could not open output file [$main_out] $!\n";
	
	foreach my $gene (keys %gene_est) {
		
		my $est = $gene_est{$gene};
		
		print MAINOUT"$gene\t$est\n";
		
	}
	
	close MAINOUT;
}

exit;
