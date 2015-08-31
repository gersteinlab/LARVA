##############################################################################
Large-scale Analysis of Recurrent Variants in Annotations (LARVA)
CODE README FILE

LARVA version 2.0
July 15, 2015

Lucas Lochovsky	and Jing Zhang
Gerstein Lab
Yale University
##############################################################################

Contents:
A) Prerequisite Software
B) File List
C) Build Instructions
D) Prerequisite User-supplied Data
E) LARVA Data Context
F) Usage of "larva" (master program)

(A) Prerequisite Software

LARVA is developed on Unix-based operating systems like Linux and Mac OS X. These types of operating systems are fully supported. It is possible that LARVA may function under Windows in an environment like Cygwin, but this has not been tested.

The following software are required to run LARVA. The "version tested" fields indicate the versions of each dependency that have been tested with LARVA. Earlier versions may work, but are unsupported.

1) gcc: The GNU project's C and C++ compiler. Necessary to produce executables from the C++ source code in the LARVA distribution.
	Link: http://gcc.gnu.org/
	Version tested: 4.4.6
	
2) make: The GNU make utility, used to automate compilation of C++ source code in LARVA.
	Link: https://www.gnu.org/software/make/
	Version tested: 3.81
	
3) bigWigAverageOverBed: Required to compute the DNA replication timing signal of each genome region under study for recurrent mutation. Must be in the same directory as LARVA or it will fail.
	*** THE 64-BIT LINUX VERSION OF THIS SCRIPT IS INCLUDED WITH LARVA. FOR OTHER VERSIONS USE THE FOLLOWING LINK.
	Link: http://genome.ucsc.edu/goldenpath/help/bigWig.html (scroll to end of page)
	Version tested: 2.0
	
4) BEDTools: Required for computing intersections between annotations and blacklist regions.
	Link: https://github.com/arq5x/bedtools2
	Version tested: 2.13.3
	
(B) File list

1) annotations: The folder for the prerequisite data context.

2) bigWigAverageOverBed: A utility script (for 64-bit Linux only) from the UCSC Genome Browser that calculates a genome signal's (in bigWig format) average over a series of intervals (in BED format). Used by LARVA to compute the DNA replication timing of each annotation in the input annotation file.

3) dependency-check.sh: Checks that all of LARVA's prerequisite software is installed at compile time.

4) larva.cpp: The master C++ source code that computes variant intersections with annotations, and calls the other .cpp files for the significance computation.

5) makefile: The script that compiles the C++ source code.

6) moment.estimator.cpp: The C++ source code for doing model fitting.

7) moment.estimator.h: The header file to allow use of functions from "moment.estimator.cpp" in "larva.cpp"

8) p.value.calc.cpp: The C++ source code for doing p-value calculations.

9) p.value.calc.h: The header file to allow use of functions from "p.value.calc.cpp" in "larva.cpp"

10) README.txt: This file. The file to explain all other files.

11) version.h: The header file that contains the current version number of LARVA. Ensures consistent version number reporting across all LARVA files.
	
(C) Build Instructions

Before LARVA can be used, its C++ source code must be compiled into executable binaries. All the requisite commands have been collected in the "makefile" file in the LARVA "code" directory. To initiate C++ compilation, "cd" into LARVA's code directory and run the "make" command.

Command summary:

	cd [larva code directory]
	make
	
(D) Prerequisite User-supplied Data

LARVA's expects variants in an input file that contains the variants pooled from all the samples under study in interval format. Specifically, this tab-delimited format is:

(chr, start, stop, cancer type, sample name)

Additional columns after these first five are allowed. For variant files in VCF format, we recommend the use of the vcf2bed.py script for conversion, available at: https://code.google.com/p/bedops/downloads/detail?name=vcf2bed.py

(E) LARVA Data Context

The data context files must be downloaded from larva.gersteinlab.org and placed in the "code/annotations" folder to run LARVA.

To compute the DNA replication timings of each genome region for the purpose of regional mutation rate correction, a replication timing signal track in bigWig format is provided. This data was sourced from Chen et al. (PMID: 20103589), and processed into bigWig format using the utilities available for download at the following URL.
	Expected path to replication timing file: code/annotations/replication_timing.bw
	bigWig URL: http://genome.ucsc.edu/goldenpath/help/bigWig.html
	
Blacklist regions refer to genome intervals for which read mappability is extremely difficult. Many genome analyses exclude these regions from analysis due to the unreliability of results. Blacklist regions were obtained from the UCSC Genome Browser at the following URL.
	Expected path to blacklist regions file: code/annotations/blacklist_regions.bed
	Blacklist regions URL: http://genome.ucsc.edu/ --> Tables
		Group: Mapping and Sequencing
		Track: Mappability
		Table: DAC Blacklist (wgEncodeDacMapabilityConsensusExcludable)
	
Gene and pseudogene annotations were derived from the main GENCODE v15 annotation file.
	Expected path to gene annotations file: code/annotations/genes.bed
	Expected path to pseudogene annotations file: code/annotations/pgenes.bed
	GENCODE URL: http://www.gencodegenes.org/releases/15.html
	
(F) Usage of "larva" (master program)

Usage: larva -vf|--variant-file [variant file] -af|--annotation-file [annotation file] -o|--output-file [output file] [-b]

Synopsis: larva takes the variants in the [variant file] and intersects them with the annotations in the [annotation file], counting up the number of intersecting variants for each annotation. Additionally, the average replication timing is calculated for each annotation using the user provided replication timing file in the "code/annotations" folder. This program also uses gene and pseudogene annotations provided by the user to flag annotations that intersect genes and/or pseudogenes. This flag uses the following values:

0 --> no gene or pseudogene intersection
1 --> gene intersection
2 --> pseudogene intersection
3 --> both gene and pseudogene intersection

Optionally, one may run the program with the [-b] option, indicating that the program should exclude any annotations from the [annotation file] that intersect the blacklist regions file in the "code/annotations" folder.

The results are produced in the [output file] with the following tab-delimited format:

(chr, start, stop, name, mutation count, gene flag, blacklist flag, length,
DNA replication timing, BBD model p-value, binomial model p-value, BBD model with repl timing correction p-value, binomial model with repl timing correction p-value, BBD model with BH adjustment p-value, BBD model with repl timing correction and BH adjustment p-value, binomial model with BH adjustment p-value, binomial model with repl timing correction and BH adjustment p-value)

The blacklist flag is 1 if the annotation intersects a blacklist region, and 0 otherwise.
