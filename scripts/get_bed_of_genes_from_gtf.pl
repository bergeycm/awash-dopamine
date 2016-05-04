#!/usr/bin/perl

use strict;
use warnings;

use List::Util qw(min max);

# ========================================================================================
# --- Script to convert GTF refGene file to BED file
# ========================================================================================

# Read in GTF file, get all coordinates of exons 
# (put them in hash indexed by gene name containing arrays),
# and output the min and max coordinates as a BED file 

# Should not be hard-coded
my $gtf_in = shift
	or die "ERROR: Must supply input GTF file via command line.\n";
chomp $gtf_in;

open (GTF, "<$gtf_in")
	or die "ERROR: Could not open input GTF file, [$gtf_in]. $!\n";

# Create hash to hold coordinates and chromosomal locations
my %coords;
my %chrs;

while (<GTF>) {

	my @gtf_info = split /\t/, $_;
	
	my $chr   = $gtf_info[0];
	my $start = $gtf_info[3];
	my $end   = $gtf_info[4];
	my $ids   = $gtf_info[8];
	
	my @id_info = split /;/, $ids;
	my $gene_id;
	if ($id_info[0] =~ /gene_id "([^"]+)"/) {
		$gene_id = $1;
	} else {
		die "ERROR: Malformed gene_id: [" . $id_info[0] . "]\n";
	}
	
	# Store chromosome in hash
	$chrs{$gene_id} = $chr;
		
	# Put this line's coordinates in hash indexed by gene name containing arrays.
	push(@{$coords{$gene_id}}, $start);
	push(@{$coords{$gene_id}}, $end);

}

close GTF;

# For each gene, print min and max coordinate as BED file.

foreach my $gene (keys %coords) {
	
	my $min = min @{$coords{$gene}};
	my $max = max @{$coords{$gene}};

	print $chrs{$gene} . "\t" . $min . "\t" . $max . "\t" . $gene . "\n";

}

exit;
