#!/usr/bin/perl
require "rand-variant-place.pl";

my $numvar = 4;
my $qbasename = "v1";
my @names = ("A1BG", "A1CF", "A2M", "A2ML1", "A4GALT");
my @cumul_prob = (41707108, 162783527, 484016454, 748202197, 788983933);
my $area = 788983933;
my %regions;
$regions{"A1BG"} = "chr19\t58864770\t58864803";
$regions{"A1CF"} = "chr10\t52619602\t52619700";
$regions{"A2M"} = "chr12\t9243797\t9244025";
$regions{"A2ML1"} = "chr12\t8975248\t8975309";
$regions{"A4GALT"} = "chr22\t43088899\t43089957";

my $names_ref = \@names;
my $cumul_prob_ref = \@cumul_prob;
my $regions_ref = \%regions;

my $return_array_ref = rand_variant_place($numvar, $qbasename, $names_ref, $cumul_prob_ref, $regions_ref, $area);
exit();
