#!/usr/bin/perl
require "rand-variant-place.pl";

# This subroutine is a wrapper for "rand-variant-place.pl" which produces multiple
# random variant files. It is designed to work with the parallel LARVA rand-pipeline

# The annotation file and database info are to be hardcoded based upon the $aropt
# option, which will determine all of these.

sub rand_variant_place_multiple {
	
	# A file containing the list of variant files. We just need the number of variants
	# in each file, and where they are in the filesystem, for creating the random files
	my $vfiles = $_[0];

	# Start index for rand generation
	# my $start = $_[1];

	# End index for rand generation
	# my $finish = $_[2];
	
	# A reference to the array with the annotation names
	my $names_ref = $_[3];
	
	# A reference to the array with the combined cumulative probability distribution
	my $cumul_prob_ref = $_[4];
	
	# A reference to the hash with the regions
	my $regions_ref = $_[5];
	
	# The total area under the curve of the cumulative probability distribution
	my $area = $_[6];
	
	# Return values
	# An array with the addresses of the @var_array and the %var_array_samps
	my @return_array;
	my @var_array = ();
	my %var_array_samps;
	$return_array[0] = \@var_array;
	$return_array[1] = \%var_array_samps;

	open VFILES, "<$vfiles" or die "Can't open $vfiles: $!\n";

	while (my $qline = <VFILES>) {
		chomp($qline);
	
		# Extract the basename from $qline
		my $qbasename;
		if (index($qline, "/") == -1) {
			$qbasename = $qline;
		} else {
			$qbasename = substr($qline, rindex($qline, "/")+1);
		}
	
		# Remove the file extension, if there is one
		# Is there a period at the 4th or 5th last position?
		if (substr($qbasename, scalar($qbasename)-4, 1) eq ".") {
			# Then remove the 3-letter file extension
			$qbasename = substr($qbasename, 0, scalar($qbasename)-4);
		} elsif (substr($qbasename, scalar($qbasename)-5, 1) eq ".") {
			# Then remove the 4-letter file extension
			$qbasename = substr($qbasename, 0, scalar($qbasename)-5);
		} # Else there is no file extension to remove
	
		# Set up output directory
	# 	my $rand_dirname = $qline."_rand";
	# 	system("mkdir -p $rand_dirname");
	
		# How many variants are in this file?
		my $wc_output = `wc -l $qline`;
		$wc_output =~ m/^(\s*)(\d+)(\s+)(.+)/;
		my $numvar = $2;
		
		# Rand generation loop here
		# for (my $i = $start; $i <= $finish; $i++) {
			my $this_return_array_ref = rand_variant_place($numvar, $qbasename, $names_ref, $cumul_prob_ref, $regions_ref, $area);
			my @this_return_array = @{$this_return_array_ref};
			
			# Unpack what's in $this_return_array_ref, and add it to the @return_array
			my @this_var_array = @{$this_return_array[0]};
			my %this_var_array_samps = %{$this_return_array[1]};
			
			# DEBUG: check the values of @this_var_array and %this_var_array_samps
# 			print "DEBUG: this_var_array: @this_var_array\n";
# 			print "DEBUG: this_var_array_samps\n";
# 			while (my ($key, $value) = each(%this_var_array_samps)) {
# 				print "$key => $value\n";
# 			}
			
			foreach $var_ref (@this_var_array) {
				my $this_chr = ${$var_ref}[0];
				my $this_start = ${$var_ref}[1];
				my $this_end = ${$var_ref}[2];
				
				my @this_var = ($this_chr, $this_start, $this_end);
				push(@var_array, \@this_var);
			}
			
			while (my ($key, $value) = each(%this_var_array_samps)) {
				if (!(exists($var_array_samps{$key}))) {
					$var_array_samps{$key} = $value;
				} else {
					$var_array_samps{$key} .= "\t".$value;
				}
			}
			# DEBUG: check the output of the var_array_samps
# 			print "DEBUG var_array_samps:\n";
# 			while (my ($key, $value) = each(%var_array_samps)) {
# 				print "$key => $value\n";
# 			}
		# }
	}
	return \@return_array;
}
1;
