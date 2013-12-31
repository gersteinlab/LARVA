#!/usr/bin/perl

# For preprocessing whole genome annotation sets for LARVA-SAM, this script
# takes the DNA replication timing data from Chen et al., remapped to hg19
# genomic intervals (first input). Taking a whole genome annotation set for
# LARVA study as a second input, this script will output, for each annotation
# from input 2, a reptiming value based on the data in input 1

# The DNA replication timing data. Assume sorted by sortBed.
# Format: (chr, start, end, reptiming)
my $reptiming_file;

# A whole genome annotation set file
# Format: (chr, start, end, name)
my $annotation_file;

# The output file
# Format: (chr, start, end, name, computed_reptiming)
my $outfile;

if (scalar(@ARGV) != 3) {
	print "Usage: perl dna-reptiming-preproc.pl [reptiming file] [annotation file] [outfile]\n";
	exit(1);
} else {
	$reptiming_file = shift(@ARGV);
	$annotation_file = shift(@ARGV);
	$outfile = shift(@ARGV);
}

# We're going to consider each of the $annotation_file's annotations in turn,
# and find the intersection using bedtools' intersectBed tool (assume it's installed)
open AFILE, "<$annotation_file" or die "Can't open $annotation_file: $!\n";
open OUTFILE, ">$outfile" or die "Can't open $outfile: $!\n";
while (my $line = <AFILE>) {
	chomp($line);
	my ($chr, $start, $end, $name) = split(/\t/, $line);

	# Temporary files we will create in the current wd
	my $in = "in.txt";
	my $out = "out.txt";
	
	open IN, ">$in" or die "Can't open $in: $!\n";
	print IN $line;
	close(IN);
	
	# Run the bedtools
	system("intersectBed -a $reptiming_file -b $in -wa | sort | uniq > $out");
	
	# Read in the resulting lines
	my @replines;
	open OUT, "<$out" or die "Can't open $out: $!\n";
	while (my $outline = <OUT>) {
		chomp($outline);
		push(@replines, $outline);
	}
	close(OUT);
	
	# Calculate the reptiming for this annotation: 
	# sum((# intersecting bp)*(reptiming for that region))/(total bp)
	my $computed_reptiming = 0;
	for (my $i = 0; $i < scalar(@replines); $i++) {
		my $cur_region = $replines[$i];
		my ($cur_chr, $cur_start, $cur_end, $cur_reptiming) = split(/\t/, $cur_region);
	
		# Need special handling for first and last ones
		if ($i == 0) {
			my $int_bp = $cur_end - $start + 1;
			$computed_reptiming += $int_bp*$cur_reptiming;
		} elsif ($i == scalar(@replines)-1) {
			my $int_bp = $end - $cur_start + 1;
			$computed_reptiming += $int_bp*$cur_reptiming;
		} else {
			my $int_bp = $cur_end - $cur_start + 1;
			$computed_reptiming += $int_bp*$cur_reptiming;
		}
	}
	
	# Now let's finish this off
	$computed_reptiming = $computing_reptiming/($end-$start+1);
	
	print OUTFILE $chr."\t".$start."\t".$end."\t".$name."\t".$computed_reptiming."\n";
}
close(AFILE);
close(OUTFILE);
exit();
