#!/usr/bin/perl

use strict;
use warnings;

use List::Util qw(sum min max);
use Term::ANSIColor;

my $region_size = 50000 * 2;	# Entire window (100kb)
my $step_size = 100;			# Distance between points for Lifting Over window
my $num_steps = $region_size / $step_size;

my $bedtools = $ENV{'BEDTOOLS'};
my $kent     = $ENV{'KENT'};

my $baboon_genome_fa       = $ENV{'GENOME_FA'};
my $rhesus_genome_fa       = $ENV{'RHESUS_GENOME_FA'};
my $baboon_to_rhesus_chain = $ENV{'PAP_RHE_CHAIN'};
my $rhesus_refgene         = $ENV{'RHESUS_REFGENE'};

my $possible_hits_perc_cutoff = 0.20;   # At least this many of the loci must LiftOver, 
                                        # and match to a single chromosome.
my $actual_hits_perc_cutoff   = 0.85;   # At least this many of the LiftOver'd loci
                                        # must match to a single chromosome.

my $rhesus_length_min = $region_size * 0.75;   # To be included, regions cannot 
                                               #  be less than this.
my $rhesus_length_max = $region_size * 1.25;   # To be included, regions cannot 
                                               #  be more than this.

my $in_file = shift;
chomp $in_file;

open (AVGS, "<$in_file")
	or die "ERROR: Couldn't open input file $in_file. $!\n";

# Make directory to hold results
my $analysis;
my $baboon_chr;
if ($in_file =~ /results\/([^_]+)_weighted_averages_(chr\d+)\.txt/) {
	$analysis = $1;
	$baboon_chr = $2;
} else {
	die "ERROR: Problem in non-standard name of input file [$in_file].\n";
}

# Make results directory (if it doesn't exist)
my $output_dir = "results/orthologs/$baboon_chr/";
system ("mkdir -p $output_dir");

# Figure out name of file to which to write results summary
# Delete it if it already exists
my $main_out_file = "results/orthologs/orthologs.${baboon_chr}.genes.out";
unlink $main_out_file;

my $header = <AVGS>;

# Make BED file
open (BED, ">${in_file}.bed")
	or die "ERROR: Could not make output BED file ${in_file}.bed. $!\n";

while (<AVGS>) {

	my @info = split;
	my $chr   = $info[0];
	my $pos   = $info[1];
	my $w_avg = $info[2];
	my $p_val = $info[3];
	chomp $p_val;

	# Figure out start and end of region to translate between genomes
	my $start = $pos - ($region_size / 2);
	my $end   = $pos + ($region_size / 2);

	# Process region
	print STDERR "$chr:$pos (has p-value of $p_val)\n";
		
	# Write coordinates to BED file.
	print BED "$chr\t$start\t$end\n";
		
	# Process this window
	&process_window ($chr, $pos, $start, $end, $p_val);
				
	print STDERR "=" x 80;
	print STDERR "\n";

}

close (BED);

# Clean up temporary files
print STDERR "Deleting temporary files that may still exist.\n";
system("rm -f ${in_file}*chr*.*.rheMac3.bed");
system("rm -f ${in_file}*chr*.*.rheMac3.unmapped.bed");
print STDERR "Finished clean up.\n";

exit;

sub process_window {

	# Get parameters
	my $chr      = $_[0];
	my $pos      = $_[1];
	my $start    = $_[2];
	my $end      = $_[3];
	my $p_val    = $_[4];
	
	# Create BED of points spaced evenly throughout window to LiftOver
	open (TMP_PTS, ">${in_file}.$chr.$pos.bed")
		or die "ERROR: Could not make BED file ${in_file}.$chr.$pos.bed. $!\n";
	for (my $p = $start; $p < $end; $p += $step_size) {
		print TMP_PTS "$chr\t$p\t" . ($p + 1) . "\n";
	}
	close TMP_PTS;
	
	# LiftOver these points that cover the window
	# Syntax: liftOver oldFile map.chain newFile unMapped
	my $liftover_cmd = "${kent}/liftOver ${in_file}.$chr.$pos.bed ";
	$liftover_cmd .= "$baboon_to_rhesus_chain ${in_file}.$chr.$pos.rheMac3.bed ";
	$liftover_cmd .= "${in_file}.$chr.$pos.rheMac3.unmapped.bed > /dev/null 2>&1";
	
	system ($liftover_cmd);

	# Delete BED of points spaced evenly throughout window
	system ("rm ${in_file}.$chr.$pos.bed");
	
	# Get tally of rhesus locations by chromosome
	my $tally_cmd = "cut -f1 ${in_file}.$chr.$pos.rheMac3.bed | ";
	$tally_cmd .= "sort | uniq -c | sed 's/^\\s\\+//g'";
	my @tally = `$tally_cmd`;
	my %chr_counts;
	my $actual_hits = 0;	# Keep track of how many hits we managed to LiftOver
	
	print STDERR colored("Window [$chr:$pos]:", 'white on_blue'), "\n";
	
	foreach (@tally) {
		my ($chr_count, $this_chr) = split (/\s/, $_);
		chomp $this_chr;
		print STDERR "\t$this_chr has $chr_count matches.\n";
		$chr_counts{$this_chr} = $chr_count;
		$actual_hits += $chr_count;
	}
	
	# Abort if no matches to avoid division by zero later on.
	if ($actual_hits == 0) {
		print STDERR colored(	"No matches. Aborting to avoid division by zero later.", 
								'red'), "\n";
		return;
	}
	
	# Figure out chromosome with most matches and how many times it's matched
	keys %chr_counts;	# Reset each iterator
	my ($best_match_chr, $max_hits) = each %chr_counts;
	
	while (my ($this_chr, $this_count) = each %chr_counts) {
		if ($this_count > $max_hits) {
			$best_match_chr = $this_chr;
			$max_hits       = $this_count;
		}
	}
	
	print STDERR "\tMost frequent chromosome: $best_match_chr\n";
	print STDERR "\tMax hits: $max_hits\n";

	# Percentage of the total possible hits 
	# (What percentage of locations are on this chromosome, of all possible hits
	# including those loci that failed to LiftOver.)
	my $possible_hits_perc = $max_hits / $num_steps;
	print STDERR "\t\tPercent of all possible: $possible_hits_perc\n";

	# (What percentage of locations are on this chromosome, of all actual hits
	# that managed to LiftOver.)
	my $actual_hits_perc = $max_hits / $actual_hits;
	print STDERR "\t\tPercent of all actual: $actual_hits_perc\n";
	
	# Test to see if this window is OK... 
	my $region_OK = 1;
	
	# Can vast majority of region be LiftOver'd and does it match the same chr?
	if ($possible_hits_perc < $possible_hits_perc_cutoff) {
		$region_OK = 0;
		print STDERR colored("Region failed. Insufficient hits or ", 'red');
		print STDERR colored("insufficient hits to single chromosome.", 'red'), "\n";
	}
	
	if ($actual_hits_perc < $actual_hits_perc_cutoff) {
		$region_OK = 0;
		print STDERR colored(	"Region failed. Insufficient hits to single chromosome.", 
								'red'), "\n";
	}
	
	# Are there outliers in base position? If so, remove them.
	# Get all basepair positions
	my $bp_cmd = "grep '$best_match_chr' ${in_file}.$chr.$pos.rheMac3.bed | ";
	$bp_cmd .= "awk '{print \$3}'";
	my @base_positions = `$bp_cmd`;

	my $outlier_count = 0;

	# Find pairwise distances between points
	my $iter = 0;
	foreach my $j (@base_positions) {
		my @distances;
		foreach my $k (@base_positions) {
			my $dist = abs($j - $k);
			push @distances, $dist;
		}

		# Find median distance
		my $median_dist;
		my @distances_sort = sort {$a <=> $b} @distances;
		my $len = @distances_sort;
		if ($len % 2 == 1) {
			$median_dist = $distances_sort[int($len/2)];
		} else {
			my $sum_mids = $distances_sort[int($len/2)-1] + $distances_sort[int($len/2)];
			$median_dist = $sum_mids / 2;
		}

		# Remove if median distance to neighbor is as big as the region size 
		# (twice the maximum of what it should, given perfect match 
		# between Papio and rhesus.)
		if ($median_dist > $region_size) {
			#print STDERR "Removing outlier: Line " . ($iter + 1) . " removed ";
			#print STDERR "since median distance to neighbor is $median_dist.\n";
			
			$outlier_count++;
			
			# Remove outlier
			my $sed_cmd = "sed -i '" . ($iter + 1) . "d' ";
			$sed_cmd .= "${in_file}.$chr.$pos.rheMac3.bed";
			system ($sed_cmd);

		} else {
	    
		    # Only increase iter if we haven't deleted a line
		    $iter++;
		}
	}
	
	# Abort if region doesn't pass
	if (!$region_OK) {
	    return;
	}
	
	# If everything OK, get min and max of region in rhesus
	my $max_line = "grep '$best_match_chr' ";
	$max_line .= "${in_file}.$chr.$pos.rheMac3.bed | ";
	my $min_line = $max_line;
	
	$max_line .= "awk 'BEGIN {max = 0} {if (\$2>max) max=\$2} END {print max}'";

	$min_line .= "awk 'NR == 1 {line = \$0; min = \$2} ";
	$min_line .= "{if (\$2<min) min=\$2} END {print min}'";
	
	my $rhesus_window_max = `$max_line`;
	chomp $rhesus_window_max;
	my $rhesus_window_min = `$min_line`;
	chomp $rhesus_window_min;
			
	print STDERR "\tRegion in rhesus: $rhesus_window_min to $rhesus_window_max\n";
	my $rhesus_region_len = $rhesus_window_max - $rhesus_window_min;
	print STDERR "\t\tRhesus region length: $rhesus_region_len\n";

	# Abort if region length is too short or too long
	if ($rhesus_region_len < $rhesus_length_min) {
	    print STDERR colored("Region length too short. ", 'red');
	    print STDERR colored("Num. of outliers: $outlier_count. Skipping.", 'red'), "\n";
	    $region_OK = 0;
	} elsif ($rhesus_region_len > $rhesus_length_max) {
	    print STDERR colored("Region length too long. Skipping.", 'red'), "\n";
	    $region_OK = 0;
	}

	# Delete LiftOver output files
	system ("rm ${in_file}.$chr.$pos.rheMac3.bed");
	system ("rm ${in_file}.$chr.$pos.rheMac3.unmapped.bed");
	
	# Abort if region doesn't pass (due to length)
	if (!$region_OK) {
            return;
        }

	# Write region to temporary BED file
	my $tmp_bed = "${in_file}.$chr.$pos.rheMac3.bed.tmp";
	
	open (TMP_RHESUS, ">$tmp_bed")
		or die "ERROR: Could not make temporary BED file $tmp_bed. $!\n";
	print TMP_RHESUS $best_match_chr . "\t";
	print TMP_RHESUS $rhesus_window_min . "\t" . $rhesus_window_max . "\n";
	close TMP_RHESUS;
	
	# Get genes (from refGene) in region
	my $results_bed = $output_dir . "$chr.$pos.rheMac3.genes.bed";
	my $bedtools_cmd = "${bedtools}/intersectBed -wa -a $rhesus_refgene ";
	$bedtools_cmd .= "-b $tmp_bed > $results_bed";
	
	system ($bedtools_cmd);
	
	# Delete temporary BED file containing region in rhesus
	system ("rm $tmp_bed");
	
	# Grab gene names from output
	my $grep_cmd = "grep -oP 'gene_name \"[A-z0-9]+\"' $results_bed | ";
	$grep_cmd .= "sort | uniq | cut -d' ' -f2";
	
	my @genes = `$grep_cmd`;
	
	# Also write it to permanent record linking the region in Papio and rhesus
	open (MAIN_OUT, ">>$main_out_file")
		or die "ERROR: Couldn't open main output file $main_out_file. $!\n";

	# Print baboon region info
	print MAIN_OUT $chr . ":";
	print MAIN_OUT $start . "-" . $end . "\t";
	
	# Print rhesus region info
	print MAIN_OUT $best_match_chr . ":";
	print MAIN_OUT $rhesus_window_min . "-" . $rhesus_window_max . "\t";
	
	# Skip printing p-value, so that we can use this output file for all statistics
	###		# Print p-value
	###		print MAIN_OUT $p_val . "\t";
	
	# Print genes in region
	foreach (@genes) {
		s/[\n"]//g;
		print MAIN_OUT $_ . ";";
	}
	print MAIN_OUT "\n";

	close MAIN_OUT;

}
