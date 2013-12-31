#!/usr/bin/perl
use DBI;
use strict;

# This is the LARVA-AIM script that works with pathway and network output from
# LARVA-Core to answer additional questions about the mutation patterns of those
# pathways/networks

# The sqlite database produced as output from a LARVA-Core run using pathways
# or networks
my $db;

# Do you want to do an afile analysis (-af):
# The afile data is dominant, with added data on how many recurrently mutated
# annotations are in each afile
# Or do you want to do an annotation analysis (-an):
# The annotation data is dominant, with added data on how many afiles each
# annotation appears in
my $opt;

# The output file
my $outfile;

if (scalar(@ARGV) != 3) {
	print "Usage: perl larva-aim.pl [database] [analysis option] [outfile]\n";
	exit(1);
} else {
	$db = shift(@ARGV);
	$opt = shift(@ARGV);
	$outfile = shift(@ARGV);
}

# Connect to $db
my $dbh = DBI->connect("dbi:SQLite:dbname=$db",'','') or die "Connection Error: $DBI::errstr\n";

if ($opt eq "-af") {

	$dbh->do("DROP TABLE IF EXISTS afile_aim");
	$dbh->do("CREATE TABLE afile_aim (afile_name varchar(50), count int)");
	$dbh->do("INSERT INTO afile_aim SELECT afile_name, count(1) FROM annotation_summary WHERE ann_nsamp > 1 GROUP BY afile_name");
	
	$dbh->do("DROP TABLE IF EXISTS afile_aim_2");
	$dbh->do("CREATE TABLE afile_aim_2 (afile_name varchar(50), nsamp int, nannot int, nvar int, count int)");
	$dbh->do("SELECT t1.*, t2.count FROM afile_summary AS t1 JOIN afile_aim AS t2 ON t1.afile_name=t2.afile_name ORDER BY t1.nsamp, t1.nannot, t1.nvar");
	
	my $sth = $dbh->prepare("SELECT * FROM afile_aim_2");
	$sth->execute or die "SQL Error: $DBI::errstr\n";
	open OUTFILE, ">$outfile" or die "Can't open $outfile: $!\n";
	
	# First, print the header row
	print OUTFILE "afile"."\t"."nsamp"."\t"."nannot"."\t"."nvar"."\t"."num_recurrent_annotations"."\n";
	
	while (my @result_row = $sth->fetchrow_array) {
		for (my $j = 0; $j < scalar(@result_row); $j++) {
			if ($j == scalar(@result_row) - 1) {
				print OUTFILE "$result_row[$j]\n";
			} else {
				print OUTFILE "$result_row[$j]\t";
			}
		}
	}
	
} elsif ($opt eq "-an") {
	
	$dbh->do("DROP TABLE IF EXISTS ann_aim");
	$dbh->do("CREATE TABLE ann_aim (chr varchar(12), start int, end int, ann_name varchar(50), nsamp int, nvar int, count int)");
	$dbh->do("SELECT chr, start, end, ann_name, nsamp, nvar, count(1) FROM annotation_summary GROUP BY chr, start, end, ann_name, nsamp, nvar ORDER BY ann_nsamp, ann_nvar");
	
	my $sth = $dbh->prepare("SELECT * FROM ann_aim");
	$sth->execute or die "SQL Error: $DBI::errstr\n";
	open OUTFILE, ">$outfile" or die "Can't open $outfile: $!\n";
	
	# First, print the header row
	print OUTFILE "chr"."\t"."start"."\t"."end"."\t"."name"."\t"."nsamp"."\t"."nvar"."\t"."num_afiles"."\n";
	
	while (my @result_row = $sth->fetchrow_array) {
		for (my $j = 0; $j < scalar(@result_row); $j++) {
			if ($j == scalar(@result_row) - 1) {
				print OUTFILE "$result_row[$j]\n";
			} else {
				print OUTFILE "$result_row[$j]\t";
			}
		}
	}
	
} else { # This is not a valid option
	print "Invalid analysis option: $opt. Must be either \'-af\' or \'-an\'.\n";
	exit(1);
}

# Verdun!
close(OUTFILE);
$dbh->disconnect();
exit();
