#!/usr/bin/perl
use Parallel::ForkManager;
use strict;
use Fcntl qw/ :flock /;
require "rand-variant-place-multiple.pl";
require "../larva-sub.pl";

# Master script for the generation of random sets of variants for LARVA to run
# on to determine statistical significance of results
# Redesigned to integrate properly with the redesigned, more efficient LARVA-Core

# ASSUMPTION: R is installed, and Rscript is in the environment $PATH
# ASSUMPTION: sqlite is installed, and sqlite3 is in the environment $PATH

# Usage: perl larva-sam.pl [vfiles] [afiles] [nrand] [ncpu] [aropt] [db]

# Subroutine for appending to one of the output files
# Arguments:
# $file: The file to append to
# $outstring: The string to print to the file
sub append_file {
	my $file = $_[0];
	my $outstring = $_[1];
	
	open OUTFILE, ">>$file" or die "Can't open $file: $!\n";
	flock(OUTFILE, LOCK_EX);
	print OUTFILE $outstring;
	close(OUTFILE);
}

# The file that contains the list of variant files (must be BED format)
my $vfiles;

# The file that contains the list of annotation files. This is the annotation set that
# the vfiles will be queried against.
my $afiles;

# Number of random samples to generate
my $nrand;

# Number of CPUs to use for parallel computations.
my $ncpu;

# Option for indicating what type of random dataset to make: exome (e) or 
# whole genome (g)
my $aropt;

# A directory for LARVA-SAM to write the output of the random LARVA runs
# This is a stopgap until a proper parallel interprocess communication setup
# can be implemented.
# my $out_db_dir;

# Specify the SQLite database to output the results to
my $db;

# Argument handling: you know the drill
if (scalar(@ARGV) != 6) {
	print "Usage: perl larva-sam.pl [vfiles] [afiles] [nrand] [ncpu] [aropt] [db]\n";
	exit(1);
} else {
	$vfiles = shift(@ARGV);
	$afiles = shift(@ARGV);
	$nrand = shift(@ARGV);
	$ncpu = shift(@ARGV);
	$aropt = shift(@ARGV);
	# $out_db_dir = shift(@ARGV);
	$db = shift(@ARGV);
}

# File tests
if (!(-e $vfiles)) {
	print "Vfiles file does not exist! Exiting.\n";
	exit(1);
} elsif (-s $vfiles == 0) {
	print "Vfiles file is empty! Exiting.\n";
	exit(1);
}

if (!(-e $afiles)) {
	print "Afiles file does not exist! Exiting.\n";
	exit(1);
} elsif (-s $afiles == 0) {
	print "Afiles file is empty! Exiting.\n";
	exit(1);
}

# Check we got a valid $aropt option
if ($aropt ne "e" && $aropt ne "g") {
	print "Invalid aropt option: Must be either \'e\' or \'g\'\n";
	exit(1);
}

# Get the number of CPUs for parallel runs
# if (`nproc`) {
# 	$ncpu=`nproc`;
# }

# Timing data
my $start=`date`;
print "Start datetime: $start\n";

# Before we begin, check $aropt and import the files we need for the variant placing
my $annotation_file;
my $regions_file;
if ($aropt eq "e") {
	# Set this to a temporary file to test things
	$annotation_file = "exome-prob-dist.txt";
	$regions_file = "gencode.v15.coding.cds.nodup.mutsig_int.bed";
} elsif ($aropt eq "g") {
	$annotation_file = "wg-prob-dist.txt";
	$regions_file = "None";
} else { # Invalid option
	print "Invalid option: $aropt. Use either \'e\' or \'g\'.\n";
	exit(1);
}

# Import the $annotation_file and $regions_file data into memory
# The first column will have the annotation names, the second will have the probability masses
my @names;
my @cumul_prob;

my $area = 0;

open AFILE, "<$annotation_file" or die "Can't open $annotation_file: $!\n";
while (my $line = <AFILE>) {
	chomp($line);
	
	my ($name, $prob) = split(/\t/, $line);
	push(@names, $name);
	$area += $prob;
	push(@cumul_prob, $area);
}
close(AFILE);

# Bring in the regions
# Format: (chr, start, end, name)
my %regions;
if ($aropt eq "e") {
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
} else { # Whole genome version. Just a stub, because the regions can be derived from the $annotation_file.
	$regions{"W"} = "G";
}

# Need the afile processing code from LARVA-Core here. This is so that this
# step is only done once.

# Annotation array, contains annotations of the format ($chr, $start, $end, $ann_name)
my @ann_array = ();

# A hash for mapping from annotations to their afiles
my %ann_array_afiles;

# Array of sample IDs
my @samples = ();

# Variable for main loop, not global variable
# A hash that maps afiles to the list of samples represented by intersecting variants
# Key: Afile name (abasename)
# Value: Reference to array of sample IDs
my %afile_sample_list;

# DEBUG
# print "DEBUG: Load annotation data into memory\n";

# Bring the annotation file data into memory
open AFILES, "<$afiles" or die "Can't open $afiles: $!\n";
while (my $aline = <AFILES>) {
	chomp($aline);
	
	# Extract the basename from $aline
	my $abasename;
	if (index($aline, "/") == -1) {
		$abasename = $aline;
	} else {
		$abasename = substr($aline, rindex($aline, "/")+1);
	}
	
	# Add this afile to %afile_summary
# 	my @starting_array = (0, 0, 0);
# 	$afile_summary{$abasename} = \@starting_array;
	my @starting_array_2 = ();
	$afile_sample_list{$abasename} = \@starting_array_2;
	
	open AFILE, "<$aline" or die "Can't open $aline: $!\n";
	while (my $line = <AFILE>) {
		chomp($line);
		my @annotation = split(/\t/, $line);
		
		# DEBUG
# 		print "What's in annotation?\n";
# 		foreach $ele (@annotation) {
# 			print "$ele\n";
# 		}
		
		# Add the afile
# 		my %afile_hash;
# 		$afile_hash{$abasename} = 1;
# 		push(@annotation, \%afile_hash);
		my $key = $annotation[0]."\t".$annotation[1]."\t".$annotation[2]."\t".$annotation[3];
		if (!(exists($ann_array_afiles{$key}))) {
			$ann_array_afiles{$key} = $abasename;
		} else {
			$ann_array_afiles{$key} .= "\t".$abasename;
		}
		push(@ann_array, \@annotation);
	}
}

# Need to populate the @samples array
# Bring variant file data into memory
open VFILES, "<$vfiles" or die "Can't open $vfiles: $!\n";
while (my $vline = <VFILES>) {
	chomp($vline);
	
	# Extract the basename from $vline
	my $vbasename;
	if (index($vline, "/") == -1) {
		$vbasename = $vline;
	} else {
		$vbasename = substr($vline, rindex($vline, "/")+1);
	}
	
	# Now take the basename and derive sample name (i.e. cut off the file extension if there is one)
	# Assume 3- or 4-letter file extension at the end
	my $sample_name = $vbasename;
	
	# Is there a period at the 4th or 5th last position?
	if (substr($sample_name, scalar($sample_name)-4, 1) eq ".") {
		# Then remove the 3-letter file extension
		$sample_name = substr($sample_name, 0, scalar($sample_name)-4);
	} elsif (substr($sample_name, scalar($sample_name)-5, 1) eq ".") {
		# Then remove the 4-letter file extension
		$sample_name = substr($sample_name, 0, scalar($sample_name)-5);
	} # Else there is no file extension to remove
	
	push(@samples, $sample_name);
}

# Divide $nrand across $ncpu
my $stepsize;
my $max;
if ($nrand >= $ncpu) {
	$stepsize = int($nrand/$ncpu);
	$max = $ncpu;
} else {
	$stepsize = 1;
	$max = $nrand;
}
my $pm = new Parallel::ForkManager($ncpu);

# system("mkdir -p $out_db_dir");

# Set up the temporary files
my $afile_summary_file = "afile_summary_file.txt";
my $annotation_summary_file = "annotation_summary_file.txt";
# my $afile_summary_file = "/nfs/storage15/ll426-data/rand-home/afile_summary_file.txt";
# my $annotation_summary_file = "/nfs/storage15/ll426-data/rand-home/annotation_summary_file.txt";

# Scrub if necessary
if (-e $afile_summary_file) {
	system("rm $afile_summary_file");
}
if (-e $annotation_summary_file) {
	system("rm $annotation_summary_file");
}

# DEBUG
# print "Parallel random variant generation and LARVA begun\n";
for (my $i = 0; $i < $max; $i++) {
	$pm->start and next; # do the fork, because chopsticks are so difficult

	my $start = ($stepsize*$i)+1;
	my $finish;
	if ($i == ($max-1)) { # Last worker
		$finish = $nrand;
	} else {
		$finish = $stepsize*($i+1);
	}
	
	# DEBUG: check that $start and $finish are correct
# 	print "$i: $start\n";
# 	print "$i: $finish\n";
	
	for (my $j = $start; $j <= $finish; $j++) {
		my $ref_array_ref = rand_variant_place_multiple($vfiles, $start, $finish, \@names, \@cumul_prob, \%regions, $area);
		my @ref_array = @{$ref_array_ref};
		my @var_array = @{$ref_array[0]};
		my %var_array_samps = %{$ref_array[1]};
		
		# DEBUG: check the output of the variant arrays
# 		print "DEBUG var_array: @var_array\n";
# 		foreach my $ele (@var_array) {
# 			print "@$ele\n";
# 		}
# 		print "DEBUG var_array_samps:\n";
# 		while (my ($key, $value) = each(%var_array_samps)) {
# 			print "$key => $value\n";
# 		}
	
		my $results_array_ref = larva(\@var_array, \%var_array_samps, \@ann_array, \%ann_array_afiles, \@samples, \%afile_sample_list);
	
		# Unpack $results_array_ref
		my @results_array = @{$results_array_ref};
		my %afile_summary = %{$results_array[0]};
		my %annotation_summary = %{$results_array[1]};
		my @mappings = @{$results_array[2]};
	
		# Set up the output db
		# my $out_db = $out_db_dir."/rand_larva_".$i;
	
		# Output the %afile_summary and the %annotation_summary hashes/tables only
		# my $dbh = DBI->connect("dbi:SQLite:dbname=$out_db",'','') or die "Connection Error: $DBI::errstr\n";
	
		# Output afile_summary file
# 		my $afile_summary_table = "afile_summary_rand_".$j;
# 		$dbh->do("DROP TABLE IF EXISTS $afile_summary_table");
# 		$dbh->do("CREATE TABLE $afile_summary_table (afile_name varchar(50), nsamp int, nannot int, nvar int)");
		while (my ($afile_name, $value) = each(%afile_summary)) {
			my @out_array = @{$value};
			my $nsamp = $out_array[0];
			my $nannot = $out_array[1];
			my $nvar = $out_array[2];
			my $outstring = $j."\t".$afile_name."\t".$nsamp."\t".$nannot."\t".$nvar."\n";
			append_file($afile_summary_file, $outstring);
			# $dbh->do("INSERT INTO $afile_summary_table VALUES (\'$afile_name\', $nsamp, $nannot, $nvar)");
		}

		# Output annotation_summary file
# 		my $annotation_summary_table = "annotation_summary_rand_".$j;
# 		$dbh->do("DROP TABLE IF EXISTS $annotation_summary_table");
# 		$dbh->do("CREATE TABLE $annotation_summary_table (chr varchar(12), start int, end int, ann_name varchar(50), nsamp int, nvar int)");
		while (my ($annotation, $value) = each(%annotation_summary)) {
			my ($chr, $start, $end, $ann_name) = split(/\t/, $annotation);
	
			my @out_array = @{$value};
			my $nsamp = $out_array[1];
			my $nvar = $out_array[2];
	
			# DEBUG
		# 	print "DEBUG: List of table outputs\n";
		# 	print "@afiles\n";
		# 	print "$chr\n";
		# 	print "$start\n";
		# 	print "$end\n";
		# 	print "$ann_name\n";
		# 	print "$nsamp\n";
		# 	print "$nvar\n";
			
			my $outstring = $j."\t".$chr."\t".$start."\t".$end."\t".$ann_name."\t".$nsamp."\t".$nvar."\n";
			append_file($annotation_summary_file, $outstring);
			# $dbh->do("INSERT INTO $annotation_summary_table VALUES (\'$chr\', $start, $end, \'$ann_name\', $nsamp, $nvar)");
		}
	}
	$pm->finish;
}
$pm->wait_all_children;

# Import all the files into the database
# First, set up the tables in the database itself using the database driver
my $dbh = DBI->connect("dbi:SQLite:dbname=$db",'','') or die "Connection Error: $DBI::errstr\n";

my $afile_summary_table = "afile_summary";
$dbh->do("DROP TABLE IF EXISTS $afile_summary_table");
$dbh->do("CREATE TABLE $afile_summary_table (rand_num int, afile_name varchar(50), nsamp int, nannot int, nvar int)");

my $annotation_summary_table = "annotation_summary";
$dbh->do("DROP TABLE IF EXISTS $annotation_summary_table");
$dbh->do("CREATE TABLE $annotation_summary_table (rand_num int, chr varchar(12), start int, end int, ann_name varchar(50), nsamp int, nvar int)");

# Prepare queries in $sqlfile
my $sqlfile = "temp.sql";
# my $sqlfile = "/nfs/storage15/ll426-data/rand-home/temp.sql";
open SQLFILE, ">$sqlfile" or die "Can't open $sqlfile: $!\n";	
print SQLFILE ".import $afile_summary_file $afile_summary_table\n";
print SQLFILE ".import $annotation_summary_file $annotation_summary_table\n";
close(SQLFILE);

# Import data
system("sqlite3 $db < $sqlfile");

# Remove temporary files
system("rm $sqlfile $afile_summary_file $annotation_summary_file");

# LARVA stats collection to be done in a separate script

# DEBUG
# print "Code's end\n";

# Timing data
my $ended=`date`;
print "End datetime: $ended\n";
exit();
