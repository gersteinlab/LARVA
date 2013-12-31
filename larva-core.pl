#!/usr/bin/perl
require "larva-sub.pl";

# Wrapper for the new larva-sub.pl subroutine version of the main LARVA code

# This script takes a list of vfiles and afiles, and performs an intersect "n"
# operation on the elements
# v2 prints everything to file, then imports into SQLite

# vfiles = variant files, a series of variant calls from samples. A single file
# with the list of files is provided, one to each line. The sample names are
# automatically derived from the base filenames.
my $vfiles;

# afiles = annotation files, a series of intervals in which we want to study
# recurrent variation patterns. A single file with the list of files is provided,
# one to each line.
my $afiles;

# The output SQLite file/database. The use of SQLite here is confined solely for
# organizing output for the user.
my $db;

# DEBUG
# my $testnum = scalar(@ARGV);
# print "$testnum\n";
# exit();

if (scalar(@ARGV) != 3) {
	print "Usage: perl larva-core.pl [vfiles] [afiles] [db]\n";
	exit(1);
} else {
	$vfiles = shift(@ARGV);
	$afiles = shift(@ARGV);
	$db = shift(@ARGV);
}

# Timing data
my $start_date = `date`;
print "Start datetime: $start_date\n";

# Data structures for starting data
# Variant array, contains variants of the format ($chr, $start, $end)
my @var_array = ();

# A hash for mapping from variants to their sample lists
my %var_array_samps;

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

# Data structures for collecting results from the walk. Written to db output at program conclusion.
# Summary stats for the afiles
# Key: afile name
# Value: Array of (nsamp, nannot, nvar)
my %afile_summary;

# Summary stats for the annotations
# Key: Annotation of format (chr, start, end, ann_name)
# Value: Array of ($nsamp, $nvar)
my %annotation_summary;

# Complete list of variant-annotation mappings
# Contains arrays of format (sampID, var_chr, var_start, var_end, afile, ann_chr, ann_start, ann_end, ann_name)
my @mappings = ();

# DEBUG
# print "DEBUG: Load variant data into memory\n";

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
	
	open VFILE, "<$vline" or die "Can't open $vline: $!\n";
	while (my $line = <VFILE>) {
		chomp($line);
		my @var = split(/\t/, $line);
		
		# DEBUG
# 		print "What's in var?\n";
# 		foreach $ele (@var) {
# 			print "$ele\n";
# 		}
		
		# Add the sample ID
		# my @samp_array = ($sample_name);
		my $key = $var[0]."\t".$var[1]."\t".$var[2];
		if (!(exists($var_array_samps{$key}))) {
			$var_array_samps{$key} = $sample_name;
		} else {
			$var_array_samps{$key} .= "\t".$sample_name;
		}
		# push(@var, \@samp_array);
		push(@var_array, \@var);
	}
}

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

my $results_array_ref = larva(\@var_array, \%var_array_samps, \@ann_array, \%ann_array_afiles, \@samples, \%afile_sample_list);

# Unpack $results_array_ref
my @results_array = @{$results_array_ref};
%afile_summary = %{$results_array[0]};
%annotation_summary = %{$results_array[1]};
@mappings = @{$results_array[2]};

# Output to database
# Open database handle
my $dbh = DBI->connect("dbi:SQLite:dbname=$db",'','') or die "Connection Error: $DBI::errstr\n";
my $sqlfile = "temp.sql";

# Three tables: afile_summary, annotation_summary, and variant-annotation mappings
# Output afile_summary table
my $afile_summary_table = "afile_summary";
$dbh->do("DROP TABLE IF EXISTS $afile_summary_table");
$dbh->do("CREATE TABLE $afile_summary_table (afile_name varchar(50), nsamp int, nannot int, nvar int)");
my $afile_summary_file = "afile_summary_file.txt";
open OUTFILE, ">$afile_summary_file" or die "Can't open $afile_summary_file: $!\n";
while (my ($afile_name, $value) = each(%afile_summary)) {
	my @out_array = @{$value};
	my $nsamp = $out_array[0];
	my $nannot = $out_array[1];
	my $nvar = $out_array[2];
	print OUTFILE $afile_name."\t".$nsamp."\t".$nannot."\t".$nvar."\n";
	# $dbh->do("INSERT INTO $afile_summary_table VALUES (\'$afile_name\', $nsamp, $nannot, $nvar)");
}
close(OUTFILE);

# Prepare queries in $sqlfile
open SQLFILE, ">$sqlfile" or die "Can't open $sqlfile: $!\n";	
print SQLFILE ".import $afile_summary_file $afile_summary_table\n";
close(SQLFILE);

# Import data
system("sqlite3 $db < $sqlfile");

# Output annotation_summary table
my $annotation_summary_table = "annotation_summary";
$dbh->do("DROP TABLE IF EXISTS $annotation_summary_table");
$dbh->do("CREATE TABLE $annotation_summary_table (afile_name varchar(50), chr varchar(12), start int, end int, ann_name varchar(50), nsamp int, nvar int)");
my $annotation_summary_file = "annotation_summary_file.txt";
open OUTFILE, ">$annotation_summary_file" or die "Can't open $annotation_summary_file: $!\n";
while (my ($annotation, $value) = each(%annotation_summary)) {
	my ($chr, $start, $end, $ann_name) = split(/\t/, $annotation);

	my @out_array = @{$value};
	my $nsamp = $out_array[1];
	my $nvar = $out_array[2];

	my $key = $chr."\t".$start."\t".$end."\t".$ann_name;
	my $afiles_string = $ann_array_afiles{$key};
	my @afiles = split(/\t/, $afiles_string);

	# DEBUG
# 	print "DEBUG: List of table outputs\n";
# 	print "@afiles\n";
# 	print "$chr\n";
# 	print "$start\n";
# 	print "$end\n";
# 	print "$ann_name\n";
# 	print "$nsamp\n";
# 	print "$nvar\n";

	foreach $af (@afiles) {
	# while (my ($af, $v) = each(%afiles_hash)) {
		print OUTFILE $af."\t".$chr."\t".$start."\t".$end."\t".$ann_name."\t".$nsamp."\t".$nvar."\n";
		# $dbh->do("INSERT INTO $annotation_summary_table VALUES (\'$af\', \'$chr\', $start, $end, \'$ann_name\', $nsamp, $nvar)");
	}
}
close(OUTFILE);

# Prepare queries in $sqlfile
open SQLFILE, ">$sqlfile" or die "Can't open $sqlfile: $!\n";	
print SQLFILE ".import $annotation_summary_file $annotation_summary_table\n";
close(SQLFILE);

# Import data
system("sqlite3 $db < $sqlfile");

# Output variant-annotation mappings table
my $mappings_table = "variant_annotation_mappings";
$dbh->do("DROP TABLE IF EXISTS $mappings_table");
$dbh->do("CREATE TABLE $mappings_table (sample varchar(50), var_chr varchar(12), var_start int, var_end int, afile_name varchar(50), ann_chr varchar(12), ann_start int, ann_end int, ann_name varchar(50))");
my $mappings_file = "mappings_file.txt";
open OUTFILE, ">$mappings_file" or die "Can't open $mappings_file: $!\n";
foreach $ele (@mappings) {
	my @out_array = @{$ele};
	my $sample = $out_array[0];
	my $var_chr = $out_array[1];
	my $var_start = $out_array[2];
	my $var_end = $out_array[3];
	my $afile = $out_array[4];
	my $ann_chr = $out_array[5];
	my $ann_start = $out_array[6];
	my $ann_end = $out_array[7];
	my $ann_name = $out_array[8];
	# $dbh->do("INSERT INTO $mappings_table VALUES (\'$sample\', \'$var_chr\', $var_start, $var_end, \'$afile\', \'$ann_chr\', $ann_start, $ann_end, \'$ann_name\')");
	print OUTFILE $sample."\t".$var_chr."\t".$var_start."\t".$var_end."\t".$afile."\t".$ann_chr."\t".$ann_start."\t".$ann_end."\t".$ann_name."\n";
}
close(OUTFILE);

# Prepare queries in $sqlfile
open SQLFILE, ">$sqlfile" or die "Can't open $sqlfile: $!\n";	
print SQLFILE ".import $mappings_file $mappings_table\n";
close(SQLFILE);

# Import data
system("sqlite3 $db < $sqlfile");

# Remove temporary files
system("rm $sqlfile $afile_summary_file $annotation_summary_file $mappings_file");

# Timing data
my $end_date = `date`;
print "End datetime: $end_date\n";

# Verdun!
$dbh->disconnect();
exit();
