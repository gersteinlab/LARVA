#!/usr/bin/perl
use DBI;
use strict;

# This script queries a database(s) produced by LARVA-SAM to calculate the
# statistics for select afiles and annotations

# The user puts a list of afiles and annotations to compute stats for (one per line)
# in the $list_file, then this script will take each one, find the numbers from
# the random runs and do the stats computations

# Option to select either afiles or annotations. Either "-af" or "-an"
my $opt;

# List of afiles/annotations to check
my $list_file;

# The output database from the LARVA-SAM step, containing info on random datasets
my $rand_db;

# The number of random datasets to expect
my $nrand;

# The database to retrieve the original observed data from
my $observed_db;

# The output file
my $outfile;

if (scalar(@ARGV) != 6) {
	print "Usage: perl larva-sam-query.pl [opt] [list file] [rand db] [nrand] ".
				"[observed data db] [output file]\n";
	exit(1);
} else {
	$opt = shift(@ARGV);
	$list_file = shift(@ARGV);
	$rand_db = shift(@ARGV);
	$nrand = shift(@ARGV);
	$observed_db = shift(@ARGV);
	$outfile = shift(@ARGV);
}

# Check validity of $opt
if ($opt ne "-af" && $opt ne "-an") {
	print "Invalid opt: $opt. Must be either \'-af\' or \'-an\'.\n";
	exit(1);
}

# Table name determined by whether we're looking for afiles or annotations
my $base_table_name;
my $lookup_attr;
if ($opt eq "-af") {
	$base_table_name = "afile_summary";
	$lookup_attr = "afile_name";
} elsif ($opt eq "-an") {
	$base_table_name = "annotation_summary";
	$lookup_attr = "ann_name";
}

# Get the number of databases in $rand_db
# my $db_count = `ls -l $rand_db | wc -l`;
# chomp($db_count);
# $db_count--;

# Divide $nrand across $db_count, to help with assembling the data across
# multiple databases
# my $stepsize;
# my $max;
# if ($nrand >= $db_count) {
# 	$stepsize = int($nrand/$db_count);
# 	$max = $db_count;
# } else {
# 	$stepsize = 1;
# 	$max = $nrand;
# }

# Declare $start and $finish
# my $start;
# my $finish;

# Clear the $outfile, if it exists
if (-e $outfile) {
	system("rm $outfile");
}

open LISTFILE, "<$list_file" or die "Can't open $list_file: $!\n";
while (my $line = <LISTFILE>) {
	chomp($line);
	
	# The arrays for the measures of mutation
	my @nsamp = ();
	my @nannot = ();
	my @nvar = ();
	
	# Search through the databases for the answers
	# Database index
	# my $i = 0;
	
	# Initialize $start and $finish
# 	$start = ($stepsize*$i)+1;
# 	if ($i == ($max-1)) { # Last worker
# 		$finish = $nrand;
# 	} else {
# 		$finish = $stepsize*($i+1);
# 	}
	
	# For loop over the random datasets
	my $dbh = DBI->connect("dbi:SQLite:dbname=$rand_db",'','') or die "Connection Error: $DBI::errstr\n";
	for (my $j = 1; $j <= $nrand; $j++) {
		# my $db = $rand_db."/rand_larva_".$i;
		# my $table_name = $base_table_name."_rand_".$j;
		
		# DEBUG
		# print "db: $db\n";
		# print "DEBUG: query: SELECT * FROM $table_name WHERE $lookup_attr = \'$line\'\n";
		
		my $sth = $dbh->prepare("SELECT * FROM $base_table_name WHERE $lookup_attr = \'$line\'");
		$sth->execute or die "SQL Error: $DBI::errstr\n";
		my @result_row = $sth->fetchrow_array;
		
		# DEBUG
		# print "DEBUG: result_row: @result_row\n";
		
		if ($opt eq "-af") {
			my $cur_nsamp = $result_row[1];
			my $cur_nannot = $result_row[2];
			my $cur_nvar = $result_row[3];
			
			push(@nsamp, $cur_nsamp);
			push(@nannot, $cur_nannot);
			push(@nvar, $cur_nvar);
		} elsif ($opt eq "-an") {
			my $cur_nsamp = $result_row[4];
			my $cur_nvar = $result_row[5];
			
			push(@nsamp, $cur_nsamp);
			push(@nvar, $cur_nvar);
		}
		$sth->finish;
		
		# Check if we need to increment database index $i
# 		if ($j == $finish) { # Update needed
# 			$i++;
# 			$start = ($stepsize*$i)+1;
# 			if ($i == ($max-1)) { # Last worker
# 				$finish = $nrand;
# 			} else {
# 				$finish = $stepsize*($i+1);
# 			}
# 		} # Otherwise, carry on
	}
	$dbh->disconnect();
	
	# Get the original observed data
	my $obs_nsamp;
	my $obs_nannot; # Only used for afiles
	my $obs_nvar;
	
	my $dbh = DBI->connect("dbi:SQLite:dbname=$observed_db",'','') or die "Connection Error: $DBI::errstr\n";
	my $sth = $dbh->prepare("SELECT * FROM $base_table_name WHERE $lookup_attr = \'$line\'");
	$sth->execute or die "SQL Error: $DBI::errstr\n";
	my @result_row = $sth->fetchrow_array;
	
	# DEBUG
	# print "DEBUG: result_row: @result_row\n";
	
	if ($opt eq "-af") {
		$obs_nsamp = $result_row[1];
		$obs_nannot = $result_row[2];
		$obs_nvar = $result_row[3];
	} elsif ($opt eq "-an") {
		$obs_nsamp = $result_row[5];
		$obs_nvar = $result_row[6];
	}
	$sth->finish;
	$dbh->disconnect();
	
	# Do the distribution fitting and p-value derivation
	# nsamp and nvar are common to both -af and -an versions
	# Nsamp analysis
	my $randfile = "rand_nsamp.txt";
	open RANDFILE, ">$randfile" or die "Can't open $randfile: $!\n";
	my $first = 1;
	foreach my $num (@nsamp) {
		if ($first) {
			print RANDFILE "$num";
			$first = 0;
		} else {
			print RANDFILE " $num";
		}
	}
	print RANDFILE "\n";
	close(RANDFILE);
	
	# Now do the random stats analysis, and print results to a file
	my $result_string = `Rscript fit_normal_cli_v2.r $randfile $obs_nsamp`;
	open OUTFILE, ">>$outfile" or die "Can't open $outfile: $!\n";
	print OUTFILE "<-- $line NSAMP -->\n";
	print OUTFILE "$result_string";
	close(OUTFILE);
	
	# Do nannot if this is an afile list analysis
	my $randfile3 = "rand_nannot.txt";
	if ($opt eq "-af") {
		open RANDFILE, ">$randfile3" or die "Can't open $randfile3: $!\n";
		$first = 1;
		foreach my $num (@nannot) {
			if ($first) {
				print RANDFILE "$num";
				$first = 0;
			} else {
				print RANDFILE " $num";
			}
		}
		print RANDFILE "\n";
		close(RANDFILE);
		
		# Now do the random stats analysis, and print results to a file
		my $result_string_3 = `Rscript fit_normal_cli_v2.r $randfile3 $obs_nannot`;
		open OUTFILE, ">>$outfile" or die "Can't open $outfile: $!\n";
		print OUTFILE "<-- $line NANNOT -->\n";
		print OUTFILE "$result_string_3";
		close(OUTFILE);
	}
	
	# Now do nvar analysis
	my $randfile2 = "rand_nvar.txt";
	open RANDFILE, ">$randfile2" or die "Can't open $randfile2: $!\n";
	$first = 1;
	foreach my $num (@nvar) {
		if ($first) {
			print RANDFILE "$num";
			$first = 0;
		} else {
			print RANDFILE " $num";
		}
	}
	print RANDFILE "\n";
	close(RANDFILE);
	
	# Now do the random stats analysis, and print results to a file
	my $result_string_2 = `Rscript fit_normal_cli_v2.r $randfile2 $obs_nvar`;
	open OUTFILE, ">>$outfile" or die "Can't open $outfile: $!\n";
	print OUTFILE "<-- $line NVAR -->\n";
	print OUTFILE "$result_string_2";
	close(OUTFILE);
	
	# Remove temporary files
	system("rm $randfile $randfile2");
	if ($opt eq "-af") {
		system("rm $randfile3");
	}
}
close(LISTFILE);
exit();
