#!/usr/bin/perl
require "rand-variant-place.pl";

# Use the subroutine in rand-variant-place.pl to generate random variant files

# Arguments to pass to rand-variant-place subroutine
my $numvar;
my $sample;

if (scalar(@ARGV) != 2) {
	print "Usage: perl rand-variant-place-wrapper.pl [numvar] [sample]\n";
	exit(1);
} else {
	$numvar = shift(@ARGV);
	$sample = shift(@ARGV);
}

# Import the $annotation_file and $regions_file data into memory
# The first column will have the annotation names, the second will have the probability masses
my @names;
my @cumul_prob;

my $area = 0;

my $annotation_file = "exome-prob-dist.txt";

open AFILE, "<$annotation_file" or die "Can't open $annotation_file: $!\n";
while (my $line = <AFILE>) {
	chomp($line);
	
	my ($name, $prob) = split(/\t/, $line);
	push(@names, $name);
	$area += $prob;
	push(@cumul_prob, $area);
}
close(AFILE);

my $regions_file = "gencode.v15.coding.cds.nodup.mutsig_int.bed";

# Bring in the regions
# Format: (chr, start, end, name)
my %regions;
open RFILE, "<$regions_file" or die "Can't open $regions_file: $!\n";
while (my $line = <RFILE>) {
	chomp($line);
	
	my ($chr, $start, $end, $name) = split(/\t/, $line);
	my $basename = $name;
	$basename =~ s/\-\d+$//;
	if (exists $regions{$basename}) { # Add region to the end
		$regions{$basename} .= ";".$chr."\t".$start."\t".$end;
	} else { # Initialize this basename
		$regions{$basename} = $chr."\t".$start."\t".$end;
	}
}
close(RFILE);

my $return_array_ref = rand_variant_place($numvar, $sample, \@names, \@cumul_prob, \%regions, $area);
my @return_array = @{$return_array_ref};
my @var_array = @{$return_array[0]};
my %var_array_samps = %{$return_array[1]};

my $outfile = $sample.".txt";
# open OUTFILE, ">$outfile" or die "Can't open $outfile: $!\n";
foreach my $var_ref (@var_array) {
	my @var = @{$var_ref};
	# print OUTFILE $var[0]."\t".$var[1]."\t".$var[2]."\n";
	print $var[0]."\t".$var[1]."\t".$var[2]."\n";
}
# close(OUTFILE);
exit();
