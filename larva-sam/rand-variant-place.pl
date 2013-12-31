#!/usr/bin/perl
use POSIX;
use strict;

# Subroutine for doing recursive BST search through the array of annotations
# and probabilities
# Argument signature:
# $randnum: The generated random number that we want to position in the arrays
# $index: Which index in the arrays are we checking?
# $stepsize: Distance to next jump in BST search
# \@names: The pointer to the array of annotation names
# \@cumul_prob: The pointer to the array of cumulative probabilities
# Return value:
# The array index that corresponds to the annotation where $randnum falls into
sub bst_search {
	# DEBUG
	# print "DEBUG: BST search\n";

	my $randnum = $_[0];
	my $index = $_[1];
	my $stepsize = $_[2];
	my @names = @{$_[3]};
	my @cumul_prob = @{$_[4]};
	
	# Set new $stepsize
	my $new_stepsize;
	if ($stepsize == 1) {
		$new_stepsize = $stepsize;
	} else {
		$new_stepsize = floor($stepsize/2);
	}
	
	# If block that determines if we have the answer or we need a recursive run
	if ($cumul_prob[$index] < $randnum && $randnum <= $cumul_prob[$index+1]) {
		# Found the target annotation at $index+1
		return ($index+1);
	} elsif ($cumul_prob[$index] == $randnum) {
		# Found the target annotation at $index
		return $index;
	} elsif ($randnum < $cumul_prob[$index]) {
		# Must search at lower indices, with the exception if $index == 0
		if ($index == 0) {
			return $index;
		} else {
			return bst_search($randnum, ($index-$stepsize), $new_stepsize, \@names, \@cumul_prob);
		}
	} elsif ($randnum > $cumul_prob[$index+1]) {
		# Must search at higher indices
		return bst_search($randnum, ($index+$stepsize), $new_stepsize, \@names, \@cumul_prob);
	}
}

sub rand_variant_place {
	# The number of variants to produce
	my $numvar = $_[0];
	
	# The sample name
	my $sample = $_[1];
	
	# A reference to the array with the annotation names
	my $names_ref = $_[2];
	
	# A reference to the array with the combined cumulative probability distribution
	my $cumul_prob_ref = $_[3];
	
	# A reference to the hash with the regions
	my $regions_ref = $_[4];
	
	# The total area under the curve of the cumulative probability distribution
	my $area = $_[5];
	
	# Return value: a two-element array containing the references to the @var_array and %var_array_samps
	my @return_array = ();
	my @var_array = ();
	my %var_array_samps;
	$return_array[0] = \@var_array;
	$return_array[1] = \%var_array_samps;
	
	# Dereference the refs
	my @names = @{$names_ref};
	my @cumul_prob = @{$cumul_prob_ref};
	my %regions = %{$regions_ref};
	
	# DEBUG: check that the dereferenced arrays/hashes work
# 	print "DEBUG: names: @names\n";
# 	print "DEBUG: cumul_prob: @cumul_prob\n";
# 	print "DEBUG: regions:\n";
# 	while (my ($key, $value) = each(%regions)) {
# 		print "$key => $value\n";
# 	}
# 	exit();
	
	for (my $i = 0; $i < $numvar; $i++) {
		# DEBUG
		# print "DEBUG: Loop start\n";

		# Code for variant dropping
		my $randnum = int(rand($area));

		# Find which annotation this corresponds to
		# First iteration: take the array's max index div 2
		my $first_index = floor((scalar(@names)-1)/2);
		my $first_stepsize = floor((scalar(@names)-1)/4);
		
		# DEBUG
# 		print "first_index: $first_index\n";
# 		print "first_stepsize: $first_stepsize\n";
		
		my $ann_index = bst_search($randnum, $first_index, $first_stepsize, \@names, \@cumul_prob);

		# Once we have the $ann_index, we drop a variant in the annotation with uniform
		# probability over the annotation's length
		my $target_annotation = $names[$ann_index];
	# 	my $sth = $dbh->prepare("SELECT * FROM $coor_table WHERE name LIKE \"$target_annotation-%\"");
	# 	$sth->execute or die "SQL Error: $DBI::errstr\n";

		my $target_region_string = $regions{$target_annotation};
		my @target_regions = split(/;/, $target_region_string);

		# Now add up the region lengths
		my $total_length = 0;
		foreach my $target (@target_regions) {
			my ($chr, $start, $end) = split(/\t/, $target);
			$total_length += $end-$start+1;
		}

		# Pick the spot in the target regions
		my $position = int(rand($total_length));
		foreach my $target (@target_regions) {
			# DEBUG
			# print "DEBUG: Loop start\n";
		
			my ($chr, $start, $end) = split(/\t/, $target);
			my $length = $end-$start;
			if ($length < $position) {
				# DEBUG
				# print "next\n";
				
				$position -= ($length+1); # next
			} else {
				# DEBUG
				# print "last\n";
				
				my $final_pos = $start + $position;
				my $final_end = $final_pos + 1;
			
				# Output finalized variant
				my @var = ($chr, $final_pos, $final_end);
				push(@var_array, \@var);
				my $key = $chr."\t".$final_pos."\t".$final_end;
				$var_array_samps{$key} = $sample;
				last;
			}
		}
		# DEBUG: check %var_array_samps values
		# my $ref = \%var_array_samps;
# 		print "DEBUG: var_array_samps:\n";
# 		while (my ($key, $value) = each(%var_array_samps)) {
# 			print "$key => $value\n";
# 		}
	}
	
	# New approach: populate %var_array_samps after the main loop
	# foreach my $var_ref (@var_array) {
		# my @var = @{$var_ref};
		# my $key = $var[0]."\t".$var[1]."\t".$var[2];
# 		my $key = "chr110001001";
# 		$var_array_samps{"$key"} = $sample;
# 		
# 		# DEBUG: check %var_array_samps values
# 		print "DEBUG: var_array_samps:\n";
# 		print "DEBUG: key: $key\n";
# 		while (my ($key, $value) = each(%var_array_samps)) {
# 			print "$key => $value\n";
# 		}	
	# }
	
	# DEBUG: check @return_array value
# 	print "DEBUG: return_array: @return_array\n";
# 	print "DEBUG: return_array[0]: @{$return_array[0]}\n";
# 	foreach my $ele (@{$return_array[0]}) {
# 		print "@$ele\n";
# 	}
# 	print "DEBUG: return_array[1]: %{$return_array[1]}\n";
# 	while (my ($key, $value) = each(%{$return_array[1]})) {
# 		print "$key => $value\n";
# 	}
	
	return \@return_array;
}
1;
