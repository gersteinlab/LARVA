#######################################################################
# Large-scale Analysis of Recurrent Variants and Annotations (LARVA)	#
# README																															#
# 																																		#
# LARVA 1.0																														#
# Oct 29, 2013																												#
# 																																		#
# Lucas Lochovsky																											#
# PhD Candidate																												#
# Gerstein Lab																												#
# Yale University																											#
#######################################################################
build 1

Contents:
A) Description
B) Prerequisites
C) Usage of LARVA-Core (larva-core.pl)
D) Usage of LARVA-SAM (larva-sam/larva-sam.pl)
E) Usage of LARVA-AIM (larva-aim.pl)

(A) Description

LARVA is a computational framework for the identification of recurrent variants
and recurrently mutated annotations in a cohort of whole genome sequenced (WGS'ed)
disease patient samples. Recurrent variants refer to single nucleotide variants
(SNVs) that are present in multiple samples. Recurrently mutated annotations refer
to annotations that contain variants that span multiple samples. LARVA allows any
set of variant calls to be used with any genome annotation set, enabling recurrent
variation analysis of a wide range of genomic elements.

LARVA-Core performs the recurrent variant analysis, and stores its output in a
SQLite database. LARVA-SAM assesses the statistical significance of the results
of LARVA-Core by simulating the creation of variant call sets where the variant
positions are randomized. The recurrent variants and annotations of the random
call sets are compared to that of the actual call sets to ascertain statistical
significance. Finally, LARVA-AIM facilitates the study of pathway and network
recurrent variation by manipulating LARVA-Core's output to link gene and
pathway/network data in a single view.

(B) Prerequisites

The following software are required for the operation of all LARVA modules.

1) Perl
	Link: http://www.perl.org

2) DBI (Perl module): Database-independent interface for Perl.
	Link: http://search.cpan.org/~timb/DBI-1.630/DBI.pm
	
3) DBD::SQLite (Perl module): SQLite specific interface for Perl. Includes a
complete installation of SQLite, which is also a prerequisite of LARVA.
	Link: http://search.cpan.org/~ishigaki/DBD-SQLite-1.40/lib/DBD/SQLite.pm

The following additional software are required for the operation of LARVA-SAM.

1) Parallel::ForkManager (Perl module): Required to automatically distribute
parallelizable work across all available CPU cores.
	Link: http://search.cpan.org/~szabgab/Parallel-ForkManager-1.05/lib/Parallel/ForkManager.pm
	
2) R: Required for determining the statistical significance of recurrent variants.
"Rscript" must be in the environment $PATH variable. For more info on the $PATH
variable, refer to: http://www.cyberciti.biz/faq/unix-linux-adding-path/
	Link: http://www.r-project.org/
	
(C) Usage of LARVA-Core

LARVA-Core performs a recurrent variant analysis with a set of variant call files,
or vfiles, and a set of genome annotation files, or afiles. The analysis is invoked
with the "larva-core.pl" script.

Usage: perl larva-core.pl [vfiles list] [afiles list] [sqlite db]

[vfiles list]: This is the path to a file with the list of vfiles to intersect
with the afiles for recurrent variant analysis. This list file contains, one per
line, the path to each vfile. Each vfile corresponds to the variants of a single
sample, and are formatted as 3-column BED files.

[afiles list]: This is the path to a file with the list of afiles to intersect
with the vfiles for recurrent variant analysis. This list file contains, one per
line, the path to each afile. Each afile is formatted with 4 columns:
(chr, start, stop, name).

[sqlite db]: The path to the SQLite database file that LARVA-Core should use for
its output. Will be created if it doesn't already exist.

After "larva-core.pl" has finished, one may use the "larva-query.pl" script to
retrieve results from the output database. "larva-query.pl" walks the user
through the steps to retrieving data on variants, annotations, samples, or afiles.
Alternatively, one may use SQLite to directly query the database with SQL syntax.
The output database consists of three tables that summarize LARVA-Core's findings:

a) afile_summary: This table contains recurrent mutation data that pertains to
the sets of annotations contained in each afile. Its schema is:
(afile name, afile number of samples mutated (nsamp), afile number of annotations
recurrently mutated (nannot), afile number of recurrent variants (nvar))

b) annotation_summary: This table contains recurrent mutation data that pertains
to individual annotations. Its schema is:
(afile name, annotation chr, annotation start, annotation end, annotation name,
annotation number of samples mutated (nsamp), annotation number of recurrent
variants (nvar))

c) variant_annotation_mappings: This table contains pairs of intersecting variants
and annotations. Its schema is:
(sample, variant chr, variant start, variant end, afile name, annotation chr,
annotation start, annotation end, annotation name)

(D) Usage of LARVA-SAM (larva-sam/larva-sam.pl)

LARVA-SAM produces replicates of the original variant call set, in terms of number
of samples and number of variants, but randomizes the positions of the variants.
It also runs LARVA-Core on the random datasets to determine the patterns of
recurrent variation that would be observed with random variant positions. This
process is invoked with the "larva-sam.pl" script (in the larva-sam directory).

Usage: perl larva-sam.pl [vfiles list] [afiles list] [nrand] [ncpu] [aropt] [out db dir]

[vfiles list]: Identical to the parameter of the same name in "larva-core.pl".
Used to determine how many samples and variants to use in the random datasets.

[afiles list]: Identical to the parameter of the same name in "larva-core.pl".
Used for the LARVA-Core runs on the random datasets.

[nrand]: The number of random datasets to generate.

[ncpu]: The number of CPU cores LARVA-SAM should use to parallelize its workload.

[aropt]: There are two options here:
"e": Randomize variants over the human exome.
"g": Randomize variants over the whole human genome.

[out db dir]: The directory to which the SQLite databases that hold the output
are written.

After "larva-sam.pl" has finished, one may use "larva-sam-query.pl" to determine
the statistical significance of the observed recurrent variation in afiles or
annotations.

Usage: perl larva-sam-query.pl [opt] [list file] [larva db dir] [nrand]
[observed data db] [output file]

[opt]: There are two options here:
"-af": Compute statistical significance for a list of afiles
"-an": Compute statistical significance for a list of annotations

[list file]: The list of afiles or annotations one wants statistical significance
computations for. Each afile/annotation appears one per line in this file.

[larva db dir]: The directory of output databases produced by "larva-sam.pl"

[nrand]: The number of random datasets that were generated.

[observed data db]: The output database produced by "larva-core.pl" with the
actual data. Used for comparison with random data.

[output file]: Contains the output of "larva-sam-query.pl". This file contains
data on the Normal distributions fitted to each measure of mutation. Measures of
mutation here refer to:

"-af" option:
a) Afile's number of samples mutated (nsamp)
b) Afile's number of annotations recurrently mutated (nannot)
c) Afile's number of recurrent variants (nvar)

"-an" option:
a) Annotation's number of samples mutated (nsamp)
b) Annotation's number of recurrent variants (nvar)

For each afile/annotation in the [list file], and each appropriate measure of
mutation, the mean, standard deviation, the 95% confidence interval bounds
(i.e. (mean-(2*SD)) and (mean+(2*SD))) of the fitted Normal distribution is
reported. The fifth value is the p-value derived from comparing the actual, observed
data to the Normal distribution. Using this number, one can determine if an
annotation, or afile's set of annotations, are recurrently mutated significantly
more than expected.

Sample output:
(Ignore the rows that say "mean" and "sd". This is added by R and is not correct.)

Original query was for two annotations, CASR-7 and ADCY6-12.

<-- CASR-7 NSAMP -->
mean   sd mean mean      
   42.749    6.521855    29.7052890    55.7927110    0.4846501 
<-- CASR-7 NVAR -->
mean   sd mean mean      
   6.0573    2.442707    1.1718863    10.9427137    0.2132177 
<-- ADCY6-12 NSAMP -->
mean   sd mean mean      
   43.0698    6.50347    30.06285909    56.07674091    0.04436481 
<-- ADCY6-12 NVAR -->
mean   sd mean mean      
   4.541301    2.115584    0.3101326    8.7724674    0.3990282

(E) Usage of LARVA-AIM (larva-aim.pl)

LARVA-AIM is an analysis integration module that uses the LARVA-Core output of
a pathway or network analysis (where individual pathways and network pairs are
represented with afiles) to place recurrently mutated genes in their pathway or
network context.

Usage: perl larva-aim.pl [database] [outfile]

[database]: The SQLite database produced as output from a LARVA-Core run using
pathways or networks

[outfile]: The file with the output of "larva-aim.pl". The output file will
contain a tab-delimited list with the following schema:
(annotation chr, annotation start, annotation end, annotation name, annotation
number of samples mutated, annotation number of recurrent variants, number of afiles
this annotation appears in)
