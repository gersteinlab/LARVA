#!/usr/bin/perl
use DBI;

# Like the name says, this subroutine sorts intervals. This is modified from the
# original version in that we expect the two intervals $a and $b to be arrayrefs,
# rather than strings. Also, they have an arbitrary number of elements after the
# (chr, start, end).
# Useful for sorting entire arrays, or for just looking at a pair.
sub sortIntervals {
	# These two arrayrefs correspond to two intervals
	my $a = $_[0];
	my $b = $_[1];
	
	my @a_arr = @{$a};
	my @b_arr = @{$b};
	
# 	my ($a_chr, $a_start, $a_end) = split(/\t/, $a);
# 	my ($b_chr, $b_start, $b_end) = split(/\t/, $b);

	# We need to modify this because $a and $b might have an arbitrary number of elements
	my $a_chr = $a_arr[0];
	my $a_start = $a_arr[1];
	my $a_end = $a_arr[2];
	
	my $b_chr = $b_arr[0];
	my $b_start = $b_arr[1];
	my $b_end = $b_arr[2];
	
	my $a_chr_int;
	if ($a_chr eq "chrX") {
		$a_chr_int = 24;
	} elsif ($a_chr eq "chrY") {
		$a_chr_int = 25;
	} elsif ($a_chr eq "chrM") {
		$a_chr_int = 26;
	} else {
		$a_chr_int = $a_chr;
		$a_chr_int =~ s/chr//;
	}
	
	my $b_chr_int;
	if ($b_chr eq "chrX") {
		$b_chr_int = 24;
	} elsif ($b_chr eq "chrY") {
		$b_chr_int = 25;
	} elsif ($b_chr eq "chrM") {
		$b_chr_int = 26;
	} else {
		$b_chr_int = $b_chr;
		$b_chr_int =~ s/chr//;
	}
	
	if ($a_chr_int != $b_chr_int) {
		return ($a_chr_int <=> $b_chr_int);
	} elsif ($a_start != $b_start) {
		return ($a_start <=> $b_start);
	} elsif ($a_end != $b_end) {
		return ($a_end <=> $b_end);
	}
}

# Code for deduplicating the elements of an array, "borrowed" from:
# http://stackoverflow.com/questions/7651/how-do-i-remove-duplicate-items-from-an-array-in-perl
# Original name: uniq2
sub dedup {
	my %seen = ();
	my @r = ();
	foreach my $a (@_) {
		unless ($seen{$a}) {
			push @r, $a;
			$seen{$a} = 1;
		}
	}
	return @r;
}

# Code for deduplicating the elements of an array, where elements are references to intervals
sub dedup2 {
	my %seen = ();
	my @r = ();
	foreach my $a (@_) {
		my @b = @{$a};
		my $c = "";
		foreach $ele (@b) {
			$c .= "$ele\t";
		}
		$c =~ s/\t$//;
		unless ($seen{$c}) {
			push @r, $a;
			$seen{$c} = 1;
		}
	}
	return @r;
}

# This script takes a list of vfiles and afiles, and performs an intersect "n"
# operation on the elements
# v2 makes some rearrangements to what happens in the main loop.
# v3 makes some adjustments to the variant collapsing step.
# v4 redesigned watch list and intersection checking bugfix.
# v5 changed %ann_sample_list's keys to a string of the relevant annotation
# v6 changes key usage in %watch_list
# v7: another pointer overhaul, this time for the embedded @sample and %afile arrays/hashes
# sub: This version turns the code into a subroutine that expects the data to be
# passed as arguments

# vfiles = variant files, a series of variant calls from samples. A single file
# with the list of files is provided, one to each line. The sample names are
# automatically derived from the base filenames.
# my $vfiles;

# afiles = annotation files, a series of intervals in which we want to study
# recurrent variation patterns. A single file with the list of files is provided,
# one to each line.
# my $afiles;

# The output SQLite file/database. The use of SQLite here is confined solely for
# organizing output for the user.
# my $db;

# if (scalar(@ARGV) != 3) {
# 	print "Usage: perl larva-v7.pl [vfiles] [afiles] [db]\n";
# 	exit(1);
# } else {
# 	$vfiles = shift(@ARGV);
# 	$afiles = shift(@ARGV);
# 	$db = shift(@ARGV);
# }

# Argument signature: ($var_array_ref, $var_array_samps_ref, $ann_array_ref, $ann_array_afiles_ref, $samples_ref, $afile_sample_list_ref)
sub larva {

	# Data structures for starting data
	# Variant array, contains variants of the format ($chr, $start, $end)
	my $var_array_ref = $_[0];
	my @var_array = @{$var_array_ref};

	# A hash for mapping from variants to their sample lists
	my $var_array_samps_ref = $_[1];
	my %var_array_samps = %{$var_array_samps_ref};

	# Pointer to where we are in the @var_array in the walk
	my $var_pointer = 0;

	# Annotation array, contains annotations of the format ($chr, $start, $end, $ann_name)
	my $ann_array_ref = $_[2];
	my @ann_array = @{$ann_array_ref};

	# A hash for mapping from annotations to their afiles
	my $ann_array_afiles_ref = $_[3];
	my %ann_array_afiles = %{$ann_array_afiles_ref};

	# Pointer to where we are in the @ann_array in the walk
	my $ann_pointer = 0;

	# Array of sample IDs
	my $samples_ref = $_[4];
	my @samples = @{$samples_ref};

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

	# Variable for main loop, not global variable
	# A hash that maps afiles to the list of samples represented by intersecting variants
	# Key: Afile name (abasename)
	# Value: Reference to array of sample IDs
	my $afile_sample_list_ref = $_[5];
	my %afile_sample_list = %{$afile_sample_list_ref};
	
	# Need to populate the %afile_summary with the names of the afiles for keys
	while (my ($key, $value) = each(%afile_sample_list)) {
		my @starting_array = (0, 0, 0);
		$afile_summary{$key} = \@starting_array;
	}

	# DEBUG
	# print "DEBUG: Load variant data into memory\n";

	# DEBUG
	# print "DEBUG: What is the state of ann_array_afiles?\n";
	# while (my ($key, $value) = each(%ann_array_afiles)) {
	# 	print "$key => $value\n";
	# }

	# Collapse @var_array elements so that all recurrent variants (variants from different
	# samples that have overlapping coordinates) are represented by a single row

	# DEBUG
	# print "DEBUG: Variant array collapsing\n";
	# print "DEBUG: var_array before collapsing\n";
	# foreach $ele (@var_array) {
	# 	print "@$ele\n";
	# }

	# Sort all elements
	@var_array = dedup2(@var_array);
	@var_array = sort {sortIntervals($a,$b)} @var_array;

	# Sort all elements
	@ann_array = dedup2(@ann_array);
	@ann_array = sort {sortIntervals($a,$b)} @ann_array;

	# DEBUG
	# print "DEBUG: ann_array after collapsing\n";
	# foreach $ele (@ann_array) {
	# 	print "@$ele\n";
	# 	while (my ($key, $value) = each(%{${$ele}[4]})) {
	# 		print "$key => $value\n";
	# 	}
	# }

	# Begin the main loop here
	# We walk the variant and annotation lists, finding intersections, and collecting
	# information on recurrent variation for output

	# The watch list of the annotations we're currently tracking
	# Contains \@ann references in the keys
	my %watch_list;

	# A hash that ties annotations to the list of samples in which it is mutated
	# Key: Annotation reference (not the annotation array itself)
	# Value: Reference to array of sample IDs
	my %ann_sample_list;

	# The previous variant we were looking at (important to see if it was updated on next iteration)
	# Points to the index in @var_array
	my $prev_var_pointer = -1;

	# The previous annotation we were looking at (important to see if it was updated on next iteration)
	# Points to the index in @ann_array
	my $prev_ann_pointer = -1;

	# DEBUG
	# print "DEBUG: Main loop begun\n";

	# Loop termination condition: we reach the end of the variant list
	while ($var_pointer < scalar(@var_array)) {

		# DEBUG
		# print "DEBUG: Loop iteration: $var_pointer\n";
	# 	print "DEBUG: What's in watch_list?\n";
	# 	while (my ($key, $value) = each(%watch_list)) {
	# 		print "$key => $value\n";
	# 	}
	
		# Unpack current pointers
		# Variant variables
		my $cur_var_ref = $var_array[$var_pointer];
		my @cur_var = @{$cur_var_ref};
	
		my $cur_var_chr = $cur_var[0];
		my $cur_var_start = $cur_var[1];
		my $cur_var_end = $cur_var[2];
		# my $cur_var_samps_ref = $cur_var[3];
		my $key = $cur_var_chr."\t".$cur_var_start."\t".$cur_var_end;
		my $cur_var_samps_string = $var_array_samps{$key};
		my @cur_var_samps = split(/\t/, $cur_var_samps_string);
	
		# Annotation variables
		my $cur_ann_ref = $ann_array[$ann_pointer];
		my @cur_ann = @{$cur_ann_ref};
	
		# DEBUG
	# 	print "What's in cur_var?\n";
	# 	foreach $ele (@cur_var) {
	# 		print "$ele\n";
	# 	}
	# 	print "What's in cur_ann?\n";
	# 	foreach $ele (@cur_ann) {
	# 		print "$ele\n";
	# 	}
	
		my $cur_ann_chr = $cur_ann[0];
		my $cur_ann_start = $cur_ann[1];
		my $cur_ann_end = $cur_ann[2];
		my $cur_ann_name = $cur_ann[3];
		# my $cur_ann_afiles_ref = $cur_ann[4];
		my $key5 = $cur_ann_chr."\t".$cur_ann_start."\t".$cur_ann_end."\t".$cur_ann_name;
		my $cur_ann_afiles_string = $ann_array_afiles{$key5};
		# print "DEBUG: $cur_ann_afiles_string\n"; # DEBUG
		my @cur_ann_afiles = split(/\t/, $cur_ann_afiles_string);
	
		# DEBUG
	# 	if (my ($key, $value) = each(%cur_ann_afiles)) {
	# 		print "$key => $value\n";
	# 	}
	# 	print "DEBUG: What's in cur_var_samps?\n";
	# 	foreach $ele (@cur_var_samps) {
	# 		print "$ele\n";
	# 	}
	# 	print "DEBUG: What's in cur_ann_afiles?\n";
	# 	foreach $ele (@cur_ann_afiles) {
	# 		print "$ele\n";
	# 	}
	
		# Check if variant updated. If so, do watch list check
		if ($prev_var_pointer != $var_pointer) {
	
			# DEBUG
			# print "DEBUG: Variant update\n";
	
			while (my ($key, $value) = each(%watch_list)) {
	# 			my $watch_ann_ref = $key;
	# 			my @watch_ann = @{$watch_ann_ref};
	# 			
	# 			my $watch_ann_chr = $watch_ann[0];
	# 			my $watch_ann_start = $watch_ann[1];
	# 			my $watch_ann_end = $watch_ann[2];
	# 			my $watch_ann_name = $watch_ann[3];
	# 			my $watch_ann_afiles_ref = $watch_ann[4];
			
				my ($watch_ann_chr, $watch_ann_start, $watch_ann_end, $watch_ann_name) = split(/\t/, $key);
				my $key4 = $watch_ann_chr."\t".$watch_ann_start."\t".$watch_ann_end."\t".$watch_ann_name;
				my $watch_ann_afiles_string = $ann_array_afiles{$key4};
				my @watch_ann_afiles = split(/\t/, $watch_ann_afiles_string);
			
				# Keys to remove from %watch_list
				my @splice_keys = ();
			
				# Check if this new variant intersects this annotation
				if ($cur_var_chr eq $watch_ann_chr && $watch_ann_start <= $cur_var_end && $cur_var_start <= $watch_ann_end) {
			
					# DEBUG
					# print "DEBUG: Yes intersection\n";
			
					# Yes intersection: updates
				
					# Add variant-annotation pair to @mappings
					foreach $samp (@cur_var_samps) {
						foreach $afile (@watch_ann_afiles) {
						# while (my ($k, $v) = each(%cur_ann_afiles)) {
							my @push_value = ($samp, $cur_var_chr, $cur_var_start, $cur_var_end, $afile, $watch_ann_chr, $watch_ann_start, $watch_ann_end, $watch_ann_name);
							push(@mappings, \@push_value);
						}
					}
				
					# annotation.nsamp
					# Assertion: this annotation appears in %watch_list and %ann_sample_list
					# We're going to just push all the sample IDs into the array and dedup later
					foreach $samp (@cur_var_samps) {
						my $key2 = $watch_ann_chr."\t".$watch_ann_start."\t".$watch_ann_end."\t".$watch_ann_name;
						push(@{$ann_sample_list{$key2}}, $samp);
					}
				
					# annotation.nvar/afile.nvar
					# Does the @cur_var appear in multiple samples?
					if (scalar(@cur_var_samps) > 1) {
						# This annotation is already added to %annotation_summary. Increment its
						# $nvar
						my $key1 = $watch_ann_chr."\t".$watch_ann_start."\t".$watch_ann_end."\t".$watch_ann_name;
						my $sum_ref = $annotation_summary{$key1};
						${$sum_ref}[2]++;
					
						# Iterate over the keys of the $afiles hash
						foreach $afile (@watch_ann_afiles) {
						# while (my ($k, $v) = each(%watch_ann_afiles)) {
							${$afile_summary{$afile}}[2]++;
						}
					}
				
					# afile.nsamp
					foreach $samp (@cur_var_samps) {
						foreach $afile (@watch_ann_files) {
						# foreach (my ($key, $value) = each(%watch_ann_afiles)) {
					
							# DEBUG
	# 						print "DEBUG: print watch_ann_afiles\n";
	# 						print "$key => $value\n";
					
							push(@{$afile_sample_list{$afile}}, $samp);
						}
					}
							
				} else { 
					# No intersection: remove from %watch_list, and do removal scripts
					push(@splice_keys, $key);
				
					# Removal scripts (tbd)
					# For each removed annotation, do the annotation.nsamp calculation
					foreach $skey (@splice_keys) {
						# Unpack annotation
	# 					my $splice_ann_ref = $skey;
	# 					my @splice_ann = @{$splice_ann_ref};
	# 					
	# 					my $splice_ann_chr = $splice_ann[0];
	# 					my $splice_ann_start = $splice_ann[1];
	# 					my $splice_ann_end = $splice_ann[2];
	# 					my $splice_ann_name = $splice_ann[3];
	# 					my $splice_ann_afiles_ref = $splice_ann[4];
	# 					my %splice_ann_afiles = %{$splice_ann_afiles_ref};

						my ($splice_ann_chr, $splice_ann_start, $splice_ann_end, $splice_ann_name) = split(/\t/, $skey);
						my $ky = $splice_ann_chr."\t".$splice_ann_start."\t".$splice_ann_end."\t".$splice_ann_name;
						my $splice_ann_afiles_string = $ann_array_afiles{$ky};
						my @splice_ann_afiles = split(/\t/, $splice_ann_afiles_string);
					
						my $key2 = $splice_ann_chr."\t".$splice_ann_start."\t".$splice_ann_end."\t".$splice_ann_name;
						my @samp_array = @{$ann_sample_list{$key2}};
						@samp_array = dedup(@samp_array);
						${$annotation_summary{$key2}}[1] = scalar(@samp_array);
						if (scalar(@samp_array) > 1) { # This means we increment afile.nannot for each afile of the @splice_ann
							# Iterate over the %splice_ann_afiles keys
							foreach $afile (@splice_ann_afiles) {
							# while (my ($k, $v) = each(%splice_ann_afiles)) {
								${$afile_summary{$afile}}[1]++;
							}
						}
					
						# DEBUG
						# print "Deleting from watch_list\n";
					
						# At the end, actually remove element from watch list
						delete $watch_list{$skey};
						# splice(@watch_list, $sindex, 1);
					}
				}
			}
		}
		
		# We do intersection check with current variant and annotation, regardless of
		# updates. If necessary, do updates, add %watch_list, add %ann_sample_list
		# New annotations have to be added to %annotation_summary
		# Intersecting annotations are added to %watch_list and %ann_sample_list
		# if ($prev_ann_pointer != $ann_pointer) {
	
			# DEBUG
	# 		if ($prev_ann_pointer != $ann_pointer) {
	# 			print "DEBUG: Annotation update\n";
	# 		}
	
			# Must add to %annotation_summary
			my $key1 = $cur_ann_chr."\t".$cur_ann_start."\t".$cur_ann_end."\t".$cur_ann_name;
			if (!(exists($annotation_summary{$key1}))) {
				my @value1_arr = ($cur_ann_afiles_ref, 0, 0);
				$annotation_summary{$key1} = \@value1_arr;
			}
		# }
		
			# DEBUG
			# print "DEBUG: Check the values of %annotation_summary at this point\n";
	# 		while (my ($key, $value) = each(%annotation_summary)) {
	# 			print "$key => $value\n";
	# 		}
	
		# Intersection check
		if ($cur_var_chr eq $cur_ann_chr && $cur_ann_start <= $cur_var_end && $cur_var_start <= $cur_ann_end) {
			# Yes intersection
		
			# DEBUG
	# 		print "DEBUG: Yes intersection\n";
	# 		print "DEBUG: Check the cur_ann_ref value: $cur_ann_ref\n";
	# 		print "DEBUG: Check the cur_ann array: @{$cur_ann_ref}\n";
		
			# Add to %watch_list and %ann_sample_list
			my $key2 = $cur_ann_chr."\t".$cur_ann_start."\t".$cur_ann_end."\t".$cur_ann_name;
			if (!(exists($watch_list{$key2}))) {
				$watch_list{$key2} = 1;
			}
			# push(@watch_list, $cur_ann_ref);
			if (!(exists($ann_sample_list{$key2}))) {
				my @samp_array = (); # Empty array reference
				$ann_sample_list{$key2} = \@samp_array;
			}
		
			# Do updates
			# Add variant-annotation pair to @mappings
			foreach $samp (@cur_var_samps) {
				foreach $afile (@cur_ann_afiles) {
				# while (my ($k, $v) = each(%cur_ann_afiles)) {
					my @push_value = ($samp, $cur_var_chr, $cur_var_start, $cur_var_end, $afile, $cur_ann_chr, $cur_ann_start, $cur_ann_end, $cur_ann_name);
					push(@mappings, \@push_value);
				}
			}
		
			# Annotation-level updates
			# annotation.nsamp
			# We're going to just push all the sample IDs into the array and dedup later
			foreach $samp (@cur_var_samps) {
				push(@{$ann_sample_list{$key2}}, $samp);
			}
		
			# annotation.nvar/afile.nvar
			# Does the @cur_var appear in multiple samples?
			if (scalar(@cur_var_samps) > 1) {
				# This annotation is already added to %annotation_summary. Increment its
				# $nvar
				my $key1 = $cur_ann_chr."\t".$cur_ann_start."\t".$cur_ann_end."\t".$cur_ann_name;
				my $sum_ref = $annotation_summary{$key1};
				${$sum_ref}[2]++;
			
				# Iterate over the keys of the $afiles hash
				foreach $afile (@cur_ann_afiles) {
				# while (my ($k, $v) = each(%cur_ann_afiles)) {
					${$afile_summary{$afile}}[2]++;
				}
			}
		
			# afile.nsamp
			# We're going to just push all the sample IDs into the array and dedup later
			foreach $samp (@cur_var_samps) {
		
				# DEBUG
	# 			print "DEBUG: samp: $samp\n";
	# 			print "DEBUG: coor: $cur_ann_chr\t$cur_ann_start\t$cur_ann_end\n";
		
				foreach $afile (@cur_ann_afiles) {
				# foreach (my ($key, $value) = each(%cur_ann_afiles)) {
			
					# DEBUG
	# 				print "DEBUG: print cur_ann_afiles\n";
	# 				print "$afile\n";
			
					push(@{$afile_sample_list{$afile}}, $samp);
				}
			}
		} # Else no intersection
	
		# OLD CODE
	
	# Intersection check
	# If no intersection, update one list's pointers and continue
	# If intersection, then things get more interesting
	# 	if ($cur_var_chr eq $cur_ann_chr && $cur_ann_start <= $cur_var_end && $cur_var_start <= $cur_ann_end) {
	# 		# Yes intersection
	# 		
	# 		# Add to watch list, if not already present
	# 		if (!(exists($watch_list{$cur_ann_ref}))) {
	# 			push(@watch_list, $cur_ann_ref);
	# 		}
	# 		
	# 		# Iterate through watch list for intersections. Remove those for which there
	# 		# are no more intersections.
	# 		while (my ($key, $value) = each(%watch_list)) {
	# 			
	# 		
	# 	} else {
	# 		# No intersection
		
		# </OLD CODE>
		
		# Compare the two pointers and advance the pointer that is "smaller" (i.e. not
		# as far along on the walk)
		my $cmp_result = sortIntervals($cur_var_ref, $cur_ann_ref);
		if ($ann_pointer < scalar(@ann_array)-1) { # Not the last annotation
			if ($cmp_result < 0) { # Variant is smaller
				$prev_var_pointer = $var_pointer;
				$var_pointer++;
				$prev_ann_pointer = $ann_pointer;
			} else { # Annotation is smaller
				$prev_ann_pointer = $ann_pointer;
				$ann_pointer++;
				$prev_var_pointer = $var_pointer;
			}
		} else { # This is the last annotation. We only increment the $var_pointer
			$prev_var_pointer = $var_pointer;
			$var_pointer++;
			$prev_ann_pointer = $ann_pointer;
		}
		# DEBUG
	# 	print "DEBUG: What's in watch_list? (end of loop)\n";
	# 	while (my ($key, $value) = each(%watch_list)) {
	# 		print "$key => $value\n";
	# 	}
	# 	my $print_val = ${$afile_summary{"d1.txt"}}[2];
	# 	print "DEBUG: afile_summary.d1.nvar: $print_val\n";
	}

	# DEBUG
	# print "DEBUG: Main loop conclusion\n";
	# print "DEBUG: ann_sample_list contents\n";
	# while (my ($key, $value) = each(%ann_sample_list)) {
	# 	print "$key => @{$value}\n";
	# }
	# print "DEBUG: watch_list contents\n";
	# while (my ($key, $value) = each(%watch_list)) {
	# 	print "$key => $value\n";
	# }
	# print "DEBUG: afile_sample_list contents\n";
	# while (my ($key, $value) = each(%afile_sample_list)) {
	# 	print "$key => @{$value}\n";
	# }

	# Conclusion of main loop scripts
	# Wrap up anything with the last annotations in the %watch_list
	while (my ($skey, $value) = each(%watch_list)) {
		# Unpack annotation
	# 	my $splice_ann_ref = $skey;
	# 	my @splice_ann = @{$splice_ann_ref};
	# 
	# 	my $splice_ann_chr = $splice_ann[0];
	# 	my $splice_ann_start = $splice_ann[1];
	# 	my $splice_ann_end = $splice_ann[2];
	# 	my $splice_ann_name = $splice_ann[3];
	# 	my $splice_ann_afiles_ref = $splice_ann[4];
	# 	my %splice_ann_afiles = %{$splice_ann_afiles_ref};

		my ($splice_ann_chr, $splice_ann_start, $splice_ann_end, $splice_ann_name) = split(/\t/, $skey);
		my $key3 = $splice_ann_chr."\t".$splice_ann_start."\t".$splice_ann_end."\t".$splice_ann_name;
		my $splice_ann_afiles_string = $ann_array_afiles{$key3};
		my @splice_ann_afiles = split(/\t/, $splice_ann_afiles_string);
	
		# DEBUG
	# 	print "DEBUG: Check splice annotation stuff\n";
	# 	print "splice ann chr: $splice_ann_chr\n";
	# 	print "splice ann start: $splice_ann_start\n";
	# 	print "splice ann end: $splice_ann_end\n";
	# 	print "splice ann name: $splice_ann_name\n";
	
		my $key2 = $splice_ann_chr."\t".$splice_ann_start."\t".$splice_ann_end."\t".$splice_ann_name;
		# my @samp_array = @{$ann_sample_list{$key2}};
		my $samp_array_ref = $ann_sample_list{$key2};
		my @samp_array = dedup(@{$samp_array_ref});
		# push(@samp_array, "q2");
	
		# DEBUG
	# 	print "DEBUG: samp array: @{$samp_array_ref}\n";
	# 	my $int = scalar(@samp_array);
	# 	print "$int\n";
	
	# 	my $sum_ref = $annotation_summary{$key2};
	# 	${$sum_ref}[1] = scalar(@samp_array);
		${$annotation_summary{$key2}}[1] = scalar(@samp_array);
		if (scalar(@samp_array) > 1) { # This means we increment afile.nannot for each afile of the @splice_ann
			# Iterate over the %splice_ann_afiles keys
			foreach $afile (@splice_ann_afiles) {
			# while (my ($k, $v) = each(%splice_ann_afiles)) {
				${$afile_summary{$afile}}[1]++;
			}
		}
	}

	# Do the afile.nsamp calculation
	# For each $afile, deduplicate its sample list and count them
	# Iterate through %afile_sample_list keys, and update %afile_summary
	while (my ($key, $value) = each(%afile_sample_list)) {
		my @samp_array = @{$afile_sample_list{$key}};
		@samp_array = dedup(@samp_array);
		${$afile_summary{$key}}[0] = scalar(@samp_array);
	}

	# DEBUG
	# print "DEBUG: annotation_summary contents\n";
	# while (my ($key, $value) = each(%annotation_summary)) {
	# 	print "$key => @{$value}\n";
	# }
	
	# Three data structures to output:  afile_summary, annotation_summary, and variant-annotation mappings
	my @return_array = (\%afile_summary, \%annotation_summary, \@mappings);
	return \@return_array;
}
1;
