#!/usr/bin/perl
use DBI;

# This is a utility script designed to help non-technical users query an output
# database from LARVA to drill down and explore the data on recurrent variation
# encapsulated within

# This script will take user input at various stages, and walk them through the
# process of retrieving what they want.

my $query;

print "Which database would you like to query?\n";
my $db = <STDIN>;
chomp($db);

print "Would you like data from the afile/annotation summary tables (s), or data from ".
"the variant-annotation mappings table (m)?\n";
my $table = <STDIN>;
chomp($table);

# Input checking
while ($table ne "s" && $table ne "m") {
	print "That was not a valid option. Type \'s\' or\'m\'.\n";
	$table = <STDIN>;
	chomp($table);
}

if ($table eq "s") {
	print "Would you like data from the afile summary table (f), or data from ".
	"the annotation summary table (n)?\n";
	my $subtable = <STDIN>;
	chomp($subtable);
	
	# Input checking
	while ($subtable ne "f" && $subtable ne "n") {
		print "That was not a valid option. Type \'f\' or \'n\'.\n";
		$subtable = <STDIN>;
		chomp($subtable);
	}
	
	if ($subtable eq "f") {
		print "Would you like to select a single afile's data (s), or output data on all of them (a)?\n";
		my $select_opt = <STDIN>;
		chomp($select_opt);
		
		# Input checking
		while ($select_opt ne "s" && $select_opt ne "a") {
			print "That was not a valid option. Type \'s\' or \'a\'.\n";
			$select_opt = <STDIN>;
			chomp($select_opt);
		}
		
		if ($select_opt eq "s") {
		
			print "Which afile would you like?\n";
			my $afile = <STDIN>;
			chomp($afile);
		
			# Final query options: sfs [afile]
			$query = "SELECT * FROM afile_summary WHERE afile_name=\'$afile\'";
		} elsif ($select_opt eq "a") {
			# Final query options: sfa
			$query = "SELECT * FROM afile_summary";
		}
		
	} elsif ($subtable eq "n") {
		print "Would you like to select: \n".
		"(n) a single annotation, or \n".
		"(f) a single afile, or \n".
		"(a) output data on all annotations, or \n".
		"(i) specify an interval in which all contained annotations should be returned?\n";
		my $select_opt = <STDIN>;
		chomp($select_opt);
		
		# Input checking
		while ($select_opt ne "n" && $select_opt ne "f" && $select_opt ne "a" && $select_opt ne "i") {
			print "That was not a valid option. Type \'n\' or \'f\' or \'a\' or \'i\'.\n";
			$select_opt = <STDIN>;
			chomp($select_opt);
		}
		
		if ($select_opt eq "n") {
			
			# If we got into here, the user either selects an annotation by name, or
			# by its interval coordinates
			print "Please indicate if you would like to select an annotation by name (n), ".
			"or by interval coordinates (c).\n";
			my $ann_select = <STDIN>;
			chomp($ann_select);
			
			# Input checking
			while ($ann_select ne "n" && $ann_select ne "c") {
				print "That was not a valid option. Type \'n\' or \'c\'.\n";
				$ann_select = <STDIN>;
				chomp($ann_select);
			}
			
			if ($ann_select eq "n") {
				print "Which annotation would you like?\n";
				my $ann_name = <STDIN>;
				chomp($ann_name);
				
				# Final query options: snnn [ann_name]
				$query = "SELECT * FROM annotation_summary WHERE ann_name=\'$ann_name\'";
			} elsif ($ann_select eq "c") {
				print "First, specify the chromosome number\n";
				my $this_chr = <STDIN>;
				chomp($this_chr);
				
				# Input checking on chromosome
				while (!($this_chr =~ m/^chr(.+)/)) {
					print "That was not a valid chromosome. Please re-enter.\n";
					$this_chr = <STDIN>;
					chomp($this_chr);
				}
				
				print "Next, we need the start coordinate of the annotation you want.\n";
				my $this_start = <STDIN>;
				chomp($this_start);
				
				# Input checking on start coor
				while (!($this_start =~ m/(\d+)/)) {
					print "The start coordinate should be a number. Please re-enter.\n";
					$this_start = <STDIN>;
					chomp($this_start);
				}
				
				print "Finally, we need the end coordinate of the annotation you want.\n";
				my $this_end = <STDIN>;
				chomp($this_end);
				
				# Input checking on end coor
				while (!($this_end =~ m/(\d+)/)) {
					print "The end coordinate should be a number. Please re-enter.\n";
					$this_end = <STDIN>;
					chomp($this_end);
				}
				
				# Final query options: snnc [this_chr, this_start, this_end]
				$query = "SELECT * FROM annotation_summary WHERE chr=\'$this_chr\' AND start=$this_start AND end=$this_end";
			}
			
		} elsif ($select_opt eq "f") {
		
			print "Please enter the name of the afile/annotation set whose annotations you ".
			"would like to print.\n";
			my $afile = <STDIN>;
			chomp($afile);
			
			# Final query options: snf [afile]
			$query = "SELECT * FROM annotation_summary WHERE afile_name=\'$afile\'";
		
		} elsif ($select_opt eq "a") {
		
			# Final query options: sna
			$query = "SELECT * FROM annotation_summary";
		
		} elsif ($select_opt eq "i") {
			print "First, specify the chromosome number\n";
				my $this_chr = <STDIN>;
				chomp($this_chr);
				
				# Input checking on chromosome
				while (!($this_chr =~ m/^chr(.+)/)) {
					print "That was not a valid chromosome. Please re-enter.\n";
					$this_chr = <STDIN>;
					chomp($this_chr);
				}
				
				print "Next, we need the start coordinate of the range you want.\n";
				my $this_start = <STDIN>;
				chomp($this_start);
				
				# Input checking on start coor
				while (!($this_start =~ m/(\d+)/)) {
					print "The start coordinate should be a number. Please re-enter.\n";
					$this_start = <STDIN>;
					chomp($this_start);
				}
				
				print "Finally, we need the end coordinate of the range you want.\n";
				my $this_end = <STDIN>;
				chomp($this_end);
				
				# Input checking on end coor
				while (!($this_end =~ m/(\d+)/)) {
					print "The end coordinate should be a number. Please re-enter.\n";
					$this_end = <STDIN>;
					chomp($this_end);
				}
				
				# Final query options: snr [this_chr, this_start, this_end]
				$query = "SELECT * FROM annotation_summary WHERE chr=\'$this_chr\' AND start>$this_start AND end<$this_end";
			}
		}
	} elsif ($table eq "m") {
		
		print "Your options here are to: \n".
		"(a) select all annotations\n".
		"(n) select one annotation\n".
		"(v) select one variant\n".
		"(i) select all variants or annotations in a certain interval\n".
		"(s) select one sample's mappings\n".
		"(f) select one afile/annotation set's mappings\n";
		
		my $select_opt = <STDIN>;
		chomp($select_opt);
		
		# Input checking
		while ($select_opt ne "a" && $select_opt ne "n" && $select_opt ne "v" &&
					 $select_opt ne "i" && $select_opt ne "s" && $select_opt ne "f") {
			print "That was not a valid option. Type \'a\' or \'n\' or \'v\' or \'i\' or \'s\' or \'f\'.\n";
			$select_opt = <STDIN>;
			chomp($select_opt);
		}
		
		if ($select_opt eq "a") {
			# Final query options: ma
			$query = "SELECT * FROM variant_annotation_mappings";
		} elsif ($select_opt eq "n") {
			# If we got into here, the user either selects an annotation by name, or
			# by its interval coordinates
			print "Please indicate if you would like to select an annotation by name (n), ".
			"or by interval coordinates (c).\n";
			my $ann_select = <STDIN>;
			chomp($ann_select);
			
			# Input checking
			while ($ann_select ne "n" && $ann_select ne "c") {
				print "That was not a valid option. Type \'n\' or \'c\'.\n";
				$ann_select = <STDIN>;
				chomp($ann_select);
			}
			
			if ($ann_select eq "n") {
				print "Which annotation would you like?\n";
				my $ann_name = <STDIN>;
				chomp($ann_name);
				
				# Final query options: mnn [ann_name]
				$query = "SELECT * FROM variant_annotation_mappings WHERE ann_name=\'$ann_name\'";
			} elsif ($ann_select eq "c") {
				print "First, specify the chromosome number\n";
				my $this_chr = <STDIN>;
				chomp($this_chr);
				
				# Input checking on chromosome
				while (!($this_chr =~ m/^chr(.+)/)) {
					print "That was not a valid chromosome. Please re-enter.\n";
					$this_chr = <STDIN>;
					chomp($this_chr);
				}
				
				print "Next, we need the start coordinate of the annotation you want.\n";
				my $this_start = <STDIN>;
				chomp($this_start);
				
				# Input checking on start coor
				while (!($this_start =~ m/(\d+)/)) {
					print "The start coordinate should be a number. Please re-enter.\n";
					$this_start = <STDIN>;
					chomp($this_start);
				}
				
				print "Finally, we need the end coordinate of the annotation you want.\n";
				my $this_end = <STDIN>;
				chomp($this_end);
				
				# Input checking on end coor
				while (!($this_end =~ m/(\d+)/)) {
					print "The end coordinate should be a number. Please re-enter.\n";
					$this_end = <STDIN>;
					chomp($this_end);
				}
				
				# Final query options: mnc [this_chr, this_start, this_end]
				$query = "SELECT * FROM variant_annotation_mappings WHERE ann_chr=\'$this_chr\' AND ann_start=$this_start AND ann_end=$this_end";
			}
		} elsif ($select_opt eq "v") {
			print "First, specify the chromosome number\n";
			my $this_chr = <STDIN>;
			chomp($this_chr);
			
			# Input checking on chromosome
			while (!($this_chr =~ m/^chr(.+)/)) {
				print "That was not a valid chromosome. Please re-enter.\n";
				$this_chr = <STDIN>;
				chomp($this_chr);
			}
			
			print "Next, we need the start coordinate of the variant you want.\n";
			my $this_start = <STDIN>;
			chomp($this_start);
			
			# Input checking on start coor
			while (!($this_start =~ m/(\d+)/)) {
				print "The start coordinate should be a number. Please re-enter.\n";
				$this_start = <STDIN>;
				chomp($this_start);
			}
			
			print "Finally, we need the end coordinate of the variant you want.\n";
			my $this_end = <STDIN>;
			chomp($this_end);
			
			# Input checking on end coor
			while (!($this_end =~ m/(\d+)/)) {
				print "The end coordinate should be a number. Please re-enter.\n";
				$this_end = <STDIN>;
				chomp($this_end);
			}
			
			# Final query options: mv [this_chr, this_start, this_end]
			$query = "SELECT * FROM variant_annotation_mappings WHERE var_chr=\'$this_chr\' AND var_start=$this_start AND var_end=$this_end";
		} elsif ($select_opt eq "i") {
			print "First, specify the chromosome number\n";
			my $this_chr = <STDIN>;
			chomp($this_chr);
			
			# Input checking on chromosome
			while (!($this_chr =~ m/^chr(.+)/)) {
				print "That was not a valid chromosome. Please re-enter.\n";
				$this_chr = <STDIN>;
				chomp($this_chr);
			}
			
			print "Next, we need the start coordinate of the interval you want.\n";
			my $this_start = <STDIN>;
			chomp($this_start);
			
			# Input checking on start coor
			while (!($this_start =~ m/(\d+)/)) {
				print "The start coordinate should be a number. Please re-enter.\n";
				$this_start = <STDIN>;
				chomp($this_start);
			}
			
			print "Finally, we need the end coordinate of the interval you want.\n";
			my $this_end = <STDIN>;
			chomp($this_end);
			
			# Input checking on end coor
			while (!($this_end =~ m/(\d+)/)) {
				print "The end coordinate should be a number. Please re-enter.\n";
				$this_end = <STDIN>;
				chomp($this_end);
			}
			
			# Final query options: mi [this_chr, this_start, this_end]
			$query = "SELECT * FROM variant_annotation_mappings WHERE var_chr=\'$this_chr\' AND var_start>$this_start AND var_end<$this_end AND ann_chr=\'$this_chr\' AND ann_start>$this_start AND ann_end<$this_end";
		} elsif ($select_opt eq "s") {
			print "Please enter the name of the sample whose annotations you ".
			"would like to print.\n";
			my $sample = <STDIN>;
			chomp($sample);
			
			# Final query options: ms [sample]
			$query = "SELECT * FROM variant_annotation_mappings WHERE sample=\'$sample\'";
		} elsif ($select_opt eq "f") {
			print "Please enter the name of the afile/annotation set whose annotations you ".
			"would like to print.\n";
			my $afile = <STDIN>;
			chomp($afile);
			
			# Final query options: mf [afile]
			$query = "SELECT * FROM variant_annotation_mappings WHERE afile_name=\'$afile\'";
		}
	}
		
# This part comes after the query is determined. Execute the query, and produce the output.
print "One last thing: do you want the output printed to a file (f), or to the screen (s)?\n";
my $out_opt = <STDIN>;
chomp($out_opt);

while ($out_opt ne "f" && $out_opt ne "s") {
	print "That was not a valid option. Type \'f\' or \'s\'.\n";
	$out_opt = <STDIN>;
	chomp($out_opt);
}

# Open database handle
my $dbh = DBI->connect("dbi:SQLite:dbname=$db",'','') or die "Connection Error: $DBI::errstr\n";

my $sth = $dbh->prepare($query);
$sth->execute or die "SQL Error: $DBI::errstr\n";
if ($out_opt eq "f") {
	print "Name the output file\n";
	my $outfile = <STDIN>;
	chomp($outfile);
	
	open OUTFILE, ">$outfile" or die "Can't open $outfile: $!\n";
	while (my @result_row = $sth->fetchrow_array) {
		my $outstring = "";
		foreach $ele (@result_row) {
			$outstring .= "$ele\t";
		}
		$outstring =~ s/\t$//;
		$outstring .= "\n";
		print OUTFILE $outstring;
	}
	close(OUTFILE);
} elsif ($out_opt eq "s") {
	while (my @result_row = $sth->fetchrow_array) {
		my $outstring = "";
		foreach $ele (@result_row) {
			$outstring .= "$ele\t";
		}
		$outstring =~ s/\t$//;
		$outstring .= "\n";
		print $outstring;
	}
}

# Conclusion
$dbh->disconnect();
exit();
