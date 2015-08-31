#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <string>
#include <errno.h>
#include <sys/stat.h>
#include <vector>
#include <map>
#include <algorithm>
#include <time.h>
#include "version.h"
#include "p.value.calc.h"

using namespace std;

#define STRSIZE 1000

// This code performs the mutation counting step of LARVA

// INPUTS

// Variant file consisting of cancer variants pooled from all the samples/patients
// the user wishes to study

// Annotation file consisting of all the genome coordinates of a class of elements
// the user wishes to study

// Optional blacklist flag to indicate that annotations that intersect blacklist
// regions should be excluded

// DATA

// DNA replication timing bigWig file, to compute the average replication timing
// of each annotation for the purposes of local mutation rate correction
// Expected path: annotations/replication_timing.bw

// Blacklist regions file, which reflects those regions in the genome that are
// difficult to map reads to. Any annotations from the annotation file that intersect
// these are removed from the analysis and output.
// Expected path: annotations/blacklist_regions.bed

// Gene annotations file, for flagging annotation file annotations that intersect
// genes.
// Expected path: annotations/genes.bed

// Pseudogene annotations file, for flagging annotation file annotations that
// intersect pseudogenes.
// Expected path: annotations/pgenes.bed

// OUTPUTS

// Mutation counts file that lists each annotation from the input annotation file
// and includes three new columns:
// 1) the number of intersecting variants
// 2) the replication timing of this annotation
// 3) a flag that indicates whether this annotation intersects a gene or pseudogene
// 0 --> no gene or pgene intersection
// 1 --> gene intersection
// 2 --> pgene intersection
// 3 --> both intersection

int main (int argc, char* argv[]) {

	// Measure running time
	time_t t1;
	time(&t1);
	
	/* Parameters */
	// Variant file
	// Format: tab(chr, start, stop, cancer type, sample name, ...)
	string vfile;
	
	// Annotation file
	// Format: tab(chr, start, stop, name, ...)
	string afile;
	
	// Path to output file
	string outfile;
	
	/* Data files */
	
	// Location of the replication timing bigWig file
	string repfile = "annotations/replication_timing.bw";
	
	// Location of blacklist region file
	string blacklist_file = "annotations/blacklist_regions.bed";
	
	// Location of gene annotations file
	string gene_file = "annotations/genes.bed";
	
	// Location of pseudogene annotations file
	string pgene_file = "annotations/pgenes.bed";
	
	// String constants for comparisons in argument handling
	char h_string[] = "-h";
	char help_string[] = "--help";
	
	char v_string[] = "-v";
	char version_string[] = "--version";
	
	char vf_string[] = "-vf";
	char variant_file_string[] = "--variant-file";
	
	char af_string[] = "-af";
	char annotation_file_string[] = "--annotation-file";
	
	char o_string[] = "-o";
	char output_file_string[] = "--output-file";
	
 	char b_string[] = "-b";
 	int blacklist_flag = 0;
	
	// Argument handling
	if (argc == 2 && (strcmp(argv[1], h_string) == 0 || strcmp(argv[1], help_string) == 0)) {
		printf("Usage: larva -vf|--variant-file [variant file] -af|--annotation-file ");
		printf("[annotation file] -o|--output-file [output file] [-b]\n\n");
		printf("Synopsis: larva takes the variants in the [variant file] and ");
		printf("intersects them with the annotations in the [annotation file], counting up ");
		printf("the number of intersecting variants for each annotation. This is the first new ");
		printf("column added to the output. Furthermore, the user-provided data files in the \"code/annotations\" ");
		printf("folder are used to compute additional info on the [annotation file] ");
		printf("annotations. The average replication timing is calculated for each ");
		printf("annotation using the DNA replication timing file. This is ");
		printf("added to the output as a second new column. A third new column serves as ");
		printf("a flag indicating whether each annotation intersects a gene (1), a pseudogene (2) ");
		printf("neither (0), or both (3). The results are produced in the [output file]. ");
		printf("A second output file is also generated with the same name as [output file] ");
		printf("with a \"-bbd\" suffix, that contains information on the beta-binomial fit ");
		printf("and the top genome annotations with recurrent mutations. Optionally, ");
		printf("one may run the program with the [-b] option, indicating that the program ");
		printf("should exclude any annotations from the [annotation file] that intersect ");
		printf("the blacklist regions file in the \"code/annotations\" folder.\n\n");
		printf("Expected files in \"code/annotations\" folder:\n");
		printf("Replication timing file: code/annotations/replication_timing.bw\n");
		printf("Blacklist regions file: code/annotations/blacklist_regions.bed\n");
		printf("Gene annotations file: code/annotations/genes.bed\n");
		printf("Pseudogene annotations file: code/annotations/pgenes.bed\n");
		return 1;
	} else if (argc == 2 && (strcmp(argv[1], v_string) == 0 || strcmp(argv[1], version_string) == 0)) {
		printf("%s\n", vers_string);
		return 1;
	} else if (argc != 7 && argc != 8) {
		printf("Incorrect number of arguments. Use -h or --help for usage information.\n");
		return 1;
	} else {
	
		// Iterate through the argument variable array
		int counter = 1;
		while (counter < argc) {
			if (strcmp(argv[counter], vf_string) == 0 || strcmp(argv[counter], variant_file_string) == 0) {
				counter++;
				vfile = string(argv[counter]);
			} else if (strcmp(argv[counter], af_string) == 0 || strcmp(argv[counter], annotation_file_string) == 0) {
				counter++;
				afile = string(argv[counter]);
			} else if (strcmp(argv[counter], o_string) == 0 || strcmp(argv[counter], output_file_string) == 0) {
				counter++;
				outfile = string(argv[counter]);
			} else if (strcmp(argv[counter], b_string) == 0) {
				// Run the blacklist exclusion code
				blacklist_flag = 1;
			} else {
				printf("Invalid argument \'%s\'. Use -h or --help for usage information.\n", argv[counter]);
				return 1;
			}
			counter++;
		}
	}
	
	/* Check for missing arguments */
	if (vfile.empty()) {
		printf("Error: Missing vfile argument. Exiting.\n");
		return 1;
	} else if (afile.empty()) {
		printf("Error: Missing afile argument. Exiting.\n");
		return 1;
	} else if (outfile.empty()) {
		printf("Error: Missing outfile argument. Exiting.\n");
		return 1;
	}
	
	/* Input validation: file tests */
	// Variant file test
	struct stat vbuf;
	if (stat(vfile.c_str(), &vbuf)) { // Report the error and exit
		printf("Error trying to get info on variant file %s: %s\n", vfile.c_str(), strerror(errno));
		return 1;
	}
	// Check that the file is not empty
	if (vbuf.st_size == 0) {
		printf("Error: The variant file cannot be empty. Exiting.\n");
		return 1;
	}
	
	// Annotation file test
	struct stat abuf;
	if (stat(afile.c_str(), &abuf)) { // Report the error and exit
		printf("Error trying to get info on the annotation file %s: %s\n", afile.c_str(), strerror(errno));
		return 1;
	}
	// Check that the file is not empty
	if (abuf.st_size == 0) {
		printf("Error: The annotation file cannot be empty. Exiting.\n");
		return 1;
	}
	
	// Replication timing file test
	struct stat rbuf;
	if (stat(repfile.c_str(), &rbuf)) { // Report the error and exit
		printf("Error trying to get info on the replication timing file %s: %s\n", repfile.c_str(), strerror(errno));
		return 1;
	}
	// Check that the file is not empty
	if (rbuf.st_size == 0) {
		printf("Error: The replication timing file cannot be empty. Exiting.\n");
		return 1;
	}
	
	// Blacklist regions file test
	struct stat blbuf;
	if (stat(blacklist_file.c_str(), &blbuf)) { // Report the error and exit
		printf("Error trying to get info on the blacklist regions file %s: %s\n", blacklist_file.c_str(), strerror(errno));
		return 1;
	}
	// Check that the file is not empty
	if (blbuf.st_size == 0) {
		printf("Error: The blacklist regions file cannot be empty. Exiting.\n");
		return 1;
	}
	
	// Gene annotations file test
	struct stat gbuf;
	if (stat(gene_file.c_str(), &gbuf)) { // Report the error and exit
		printf("Error trying to get info on the gene annotations file %s: %s\n", gene_file.c_str(), strerror(errno));
		return 1;
	}
	// Check that the file is not empty
	if (gbuf.st_size == 0) {
		printf("Error: The gene annotations file cannot be empty. Exiting.\n");
		return 1;
	}
	
	// Pgene annotations file test
	struct stat pgbuf;
	if (stat(pgene_file.c_str(), &pgbuf)) { // Report the error and exit
		printf("Error trying to get info on the pseudogene annotations file %s: %s\n", pgene_file.c_str(), strerror(errno));
		return 1;
	}
	// Check that the file is not empty
	if (pgbuf.st_size == 0) {
		printf("Error: The pseudogene annotations file cannot be empty. Exiting.\n");
		return 1;
	}
	
	/* Check the path of the output file */
	// Extract path, if there is one
	// If there is path info in the string, there will be at least one forward slash
	unsigned int lastslash = outfile.find_last_of("/");
	if (lastslash != string::npos) { // Yes path
		string outpath = outfile.substr(0,lastslash+1);
	
		struct stat obuf;
		if (stat(outpath.c_str(), &obuf)) {
			printf("Error: The output file's directory does not exist. Exiting.\n");
			return 1;
		}
	}
	
	// Verify that bigWigAverageOverBed is in the same directory as this program
	struct stat avgbuf;
	char avgoverbed_cstr[] = "./bigWigAverageOverBed";
	if (stat(avgoverbed_cstr, &avgbuf)) {
		printf("Error: bigWigAverageOverBed is not in the same directory as \"larva\". Exiting.\n");
		return 1;
	}
	
	// Verify that BEDTools is installed
	FILE *verify_bedtools;
	char buff[STRSIZE];
	if (!(verify_bedtools = popen("command -v intersectBed", "r"))) {
		printf("Error verifying BEDTools installation. Exiting.\n");
		exit(1);
	}
	if (fgets(buff, sizeof(buff), verify_bedtools) == NULL) { // Indicates intersectBed is not available
		printf("Error: BEDTools is not installed. BEDTools must be installed to proceed. Exiting.\n");
		return 1;
	}
	pclose(verify_bedtools);
	
	/* If the blacklist flag has been used, then do intersection with blacklist
	 * regions file and only import those that are non-blacklist */
	if (blacklist_flag) {
		string bl_out = "nonblacklist.txt";
		string command = "intersectBed -a " + afile + " -b " + blacklist_file + " -v > " + bl_out;
		system(command.c_str());
		afile = bl_out;
	}
	
	/* Calculate mutation counts */
	string mut_count = "temp1.bed";
	string command2 = "intersectBed -a " + afile + " -b " + vfile + " -c > " + mut_count;
	system(command2.c_str());
	
	/* Calculate intersections with genes */
	string gene_count = "temp2.bed";
	string command3 = "intersectBed -a " + afile + " -b " + gene_file + " -c > " + gene_count;
	system(command3.c_str());
	
	/* Calculate intersections with pgenes */
	string pgene_count = "temp3.bed";
	string command4 = "intersectBed -a " + afile + " -b " + pgene_file + " -c > " + pgene_count;
	system(command4.c_str());
	
	/* Data structures */
	// The vector that contains the final output.
	// Format: (chr, start, stop, name, counts, gene flag, replication timing)
	vector<vector<string> > annotation_mutation_counts;
	
	// Bring first temp file data (mutation counts) into memory
	char linebuf_cstr[STRSIZE];
	FILE *temp1_ptr = fopen(mut_count.c_str(), "r");
	int line_num = 1;
	while (fgets(linebuf_cstr, STRSIZE-1, temp1_ptr) != NULL) {
	
		// Vector for this row
		vector<string> row;
	
		// Verify the number of columns in this file
		int numcol = 0;
		
		string linebuf = string(linebuf_cstr);
		linebuf.erase(linebuf.find_last_not_of(" \n\r\t")+1);
		unsigned int pos = linebuf.find_first_of("\t");
		while (linebuf.size() > 0) {
			
			string piece;
			if (pos != (unsigned int) string::npos) {
				piece = linebuf.substr(0,pos);
				linebuf = linebuf.substr(pos+1);
			} else {
				piece = linebuf;
				linebuf.clear();
			}
			row.push_back(piece);
			numcol++;
			
			// Format verifications
			if (numcol == 1) { // Verify chromosome format
				char c1[STRSIZE];
				if (sscanf(piece.c_str(), "chr%s", c1) != 1) {
					// if (piece.compare("chrX") != 0 && piece.compare("chrY") != 0 && piece.compare("chrM") != 0) {
						printf("Error: First column in file %s, line %d is incorrectly formatted: %s. Exiting.\n", vfile.c_str(), line_num, piece.c_str());
						return 1;
					// }
				}
			} else if (numcol == 2) { // Verify start coordinate format
				int c2;
				if (sscanf(piece.c_str(), "%d", &c2) != 1) {
					printf("Error: Second column in file %s, line %d is incorrectly formatted: %s. Exiting.\n", vfile.c_str(), line_num, piece.c_str());
					return 1;
				}
			} else if (numcol == 3) { // Verify stop coordinate format
				int c3;
				if (sscanf(piece.c_str(), "%d", &c3) != 1) {
					printf("Error: Third column in file %s, line %d is incorrectly formatted: %s. Exiting.\n", vfile.c_str(), line_num, piece.c_str());
					return 1;
				}
			}
			pos = linebuf.find_first_of("\t");
		}
		
		// Did we have too few columns?
		if (numcol < 5) {
			printf("Error: Mutation counts file %s has too few columns on line %d. Expected at least 5, but found %d. Exiting.\n", vfile.c_str(), line_num, numcol);
			return 1;
		}
		annotation_mutation_counts.push_back(row);
		line_num++;
	}
	if (!(feof(temp1_ptr)) && ferror(temp1_ptr)) { // This is an error
		char preamble[STRSIZE];
		sprintf(preamble, "There was an error reading from %s", mut_count.c_str());
		perror(preamble);
		return 1;
	}
	fclose(temp1_ptr);
	
	// Bring second temp file data (gene intersections) into memory
	FILE *temp2_ptr = fopen(gene_count.c_str(), "r");
	line_num = 1;
	while (fgets(linebuf_cstr, STRSIZE-1, temp2_ptr) != NULL) {
	
		// Verify the number of columns in this file
		int numcol = 0;
		
		string linebuf = string(linebuf_cstr);
		linebuf.erase(linebuf.find_last_not_of(" \n\r\t")+1);
		unsigned int pos = linebuf.find_first_of("\t");
		while (linebuf.size() > 0) {
		
			string piece;
			if (pos != (unsigned int) string::npos) {
				piece = linebuf.substr(0,pos);
				linebuf = linebuf.substr(pos+1);
			} else {
				piece = linebuf;
				linebuf.clear();
			}
			numcol++;
			
			if (numcol == 5) { // This is the value that indicates gene intersection
				// Convert string to integer
				int gene_flag = atoi(piece.c_str());
				if (gene_flag) { // Yes gene
					annotation_mutation_counts[line_num-1].push_back(string("1"));
				} else { // No gene
					annotation_mutation_counts[line_num-1].push_back(string("0"));
				}
			}
			
			pos = linebuf.find_first_of("\t");
		}
		
		// Did we have too few columns?
		if (numcol < 5) {
			printf("Error: Gene intersection file %s has too few columns on line %d. Expected at least 4, but found %d. Exiting.\n", afile.c_str(), line_num, numcol);
			return 1;
		}
		
		line_num++;
	}
	if (!(feof(temp2_ptr)) && ferror(temp2_ptr)) { // This is an error
		char preamble[STRSIZE];
		sprintf(preamble, "There was an error reading from %s", gene_count.c_str());
		perror(preamble);
		return 1;
	}
	fclose(temp2_ptr);
	
	// Bring third temp file data (pgene intersections) into memory
	FILE *temp3_ptr = fopen(pgene_count.c_str(), "r");
	line_num = 1;
	while (fgets(linebuf_cstr, STRSIZE-1, temp3_ptr) != NULL) {
	
		// Verify the number of columns in this file
		int numcol = 0;
		
		string linebuf = string(linebuf_cstr);
		linebuf.erase(linebuf.find_last_not_of(" \n\r\t")+1);
		unsigned int pos = linebuf.find_first_of("\t");
		while (linebuf.size() > 0) {
		
			string piece;
			if (pos != (unsigned int) string::npos) {
				piece = linebuf.substr(0,pos);
				linebuf = linebuf.substr(pos+1);
			} else {
				piece = linebuf;
				linebuf.clear();
			}
			numcol++;
			
			if (numcol == 5) { // This is the value that indicates pgene intersection
				string gene_flag_str = annotation_mutation_counts[line_num-1][5];
				int gene_flag = atoi(gene_flag_str.c_str());
				
				// Convert string to integer
				int pgene_flag = atoi(piece.c_str());
				if (pgene_flag) { // Yes pgene
					gene_flag += 2;
					char gene_flag_cstr[STRSIZE];
					sprintf(gene_flag_cstr, "%d", gene_flag);
					annotation_mutation_counts[line_num-1][5] = string(gene_flag_cstr);
				} // Else no pgene, no modifications necessary
			}
			
			pos = linebuf.find_first_of("\t");
		}
		
		// Did we have too few columns?
		if (numcol < 5) {
			printf("Error: Pgene intersection file %s has too few columns on line %d. Expected at least 3, but found %d. Exiting.\n", pgene_count.c_str(), line_num, numcol);
			return 1;
		}
		line_num++;
	}
	if (!(feof(temp3_ptr)) && ferror(temp3_ptr)) { // This is an error
		char preamble[STRSIZE];
		sprintf(preamble, "There was an error reading from %s", pgene_count.c_str());
		perror(preamble);
		return 1;
	}
	fclose(temp3_ptr);
	
	/* Replication timing computation */
	// Need to create a new BED file with unique names
	string avg_infile = "temp4.bed";
	string avg_outfile = "temp5.txt";
	int regnum = 1;
	FILE *avg_infile_ptr = fopen(avg_infile.c_str(), "w");
	for (unsigned int i = 0; i < annotation_mutation_counts.size(); i++) {
		char regnum_cstr[STRSIZE];
		sprintf(regnum_cstr, "%d", regnum);
		string outstring = annotation_mutation_counts[i][0] + "\t" + annotation_mutation_counts[i][1] + "\t" + annotation_mutation_counts[i][2] + "\t" + "reg" + string(regnum_cstr) + "\n";
		fprintf(avg_infile_ptr, outstring.c_str());
		regnum++;
	}
	fclose(avg_infile_ptr);
	
	// The actual command
	// Assumes bigWigAverageOverBed is in same directory
	string command = "./bigWigAverageOverBed " + repfile + " " + avg_infile + " " + avg_outfile;
	system(command.c_str());
	
	vector<double> reptimings;
	
	// Read the output into memory to combine with the other results
	FILE *avg_outfile_ptr = fopen(avg_outfile.c_str(), "r");
	while (fgets(linebuf_cstr, STRSIZE-1, avg_outfile_ptr) != NULL) {
		
		string linebuf = string(linebuf_cstr);
		int col_index = 0;
		while (col_index < 5) {
			unsigned int pos = linebuf.find_first_of("\t");
			linebuf = linebuf.substr(pos+1);
			col_index++;
		}
	
		// Now linebuf has the value we're looking for. Put it in the reptimings vector.
		double reptiming;
		sscanf(linebuf.c_str(), "%lf", &reptiming);
		reptimings.push_back(reptiming);
	}
	if (!(feof(avg_outfile_ptr)) && ferror(avg_outfile_ptr)) { // This is an error
		char preamble[STRSIZE];
		sprintf(preamble, "There was an error reading from %s", avg_outfile.c_str());
		perror(preamble);
		return 1;
	}
	fclose(avg_outfile_ptr);
	
	// Combine data from annotation_mutation_counts and reptimings
	for (unsigned int i = 0; i < annotation_mutation_counts.size(); i++) {
		// Convert the reptiming to a string
		char doublebuf[STRSIZE];
		sprintf(doublebuf, "%f", reptimings[i]);
		annotation_mutation_counts[i].push_back(string(doublebuf));
	}
	
	// Model fitting
	// Produce length vector
	vector<int> lengths;
	
	// Boolean for whether all lengths are the same
	bool is_unique = true;
	
	// The first observed length
	int first_length;
	
	for (unsigned int i = 0; i < annotation_mutation_counts.size(); i++) {
		int start = atoi(annotation_mutation_counts[i][1].c_str());
		int end = atoi(annotation_mutation_counts[i][2].c_str());
		int length = end - start;
		
		if (length <= 0) {
			printf("Error: Invalid length of %d in annotation file, line %d\n", length, i+1);
			printf("Length must be greater than zero\n");
			return 1;
		}
		
		if (i == 0) {
			first_length = length;
		} else {
			if (first_length != length) {
				is_unique = false;
			}
		}
		
		lengths.push_back(length);
	}
	
	// Produce counts vector
	vector<int> counts;
	
	// Boolean for whether all counts have the same value
	bool counts_unique = true;
	
	// The first observed count
	int first_count;
	
	for (unsigned int i = 0; i < annotation_mutation_counts.size(); i++) {
		int count = atoi(annotation_mutation_counts[i][4].c_str());
		
		if (i == 0) {
			first_count = count;
		} else {
			if (first_count != count) {
				counts_unique = false;
			}
		}
		
		counts.push_back(count);
	}
	
	vector<vector<double> > pval;
	vector<vector<double> > pval_corre;
	
	if (is_unique) { // Oui
		pval = p_value_equal_len(counts, lengths);
		pval_corre = pval_fixed_length_rep_time_correction(counts, lengths, reptimings);
	} else { // Non
		pval = pval_varying_length(counts, lengths);
		pval_corre = pval_varying_length_rep_time_correction(counts, lengths, reptimings);
	}
	
	vector<double> d_p_bbd = pval[1];
	vector<double> d_p_binomial = pval[0];
	vector<double> d_p_bbd_cor = pval_corre[1];
	vector<double> d_p_binomial_cor = pval_corre[0];
	
	for (unsigned int i = 0; i < d_p_bbd.size(); i++) {
		if (d_p_bbd[i] <= 0) {
			d_p_bbd[i] = 2.2e-16;
		}
	}
	
	for (unsigned int i = 0; i < d_p_bbd_cor.size(); i++) {
		if (d_p_bbd_cor[i] <= 0) {
			d_p_bbd_cor[i] = 2.2e-16;
		}
	}
	
	// Benjamini-Hochberg FDR correction
	vector<double> d_p_bbd_adj = bh_adjust(d_p_bbd);
	vector<double> d_p_binomial_adj = bh_adjust(d_p_binomial);
	vector<double> d_p_bbd_cor_adj = bh_adjust(d_p_bbd_cor);
	vector<double> d_p_binomial_cor_adj = bh_adjust(d_p_binomial_cor);
	
	// Print output
	FILE *outfile_ptr = fopen(outfile.c_str(), "w");
	
	// Print the header row first
	string header = "chr\tstart\tstop\tname\tct\tflag\trep\tlen\tp.bbd\tp.binomial\tp.bbd.cor\tp.binomial.cor\tp.bbd.adj\tp.bbd.cor.adj\tp.binomial.adj\tp.binomial.cor.adj\n";
	fprintf(outfile_ptr, header.c_str());
	
	for (unsigned int i = 0; i < annotation_mutation_counts.size(); i++) {
	
		char this_length[STRSIZE];
		sprintf(this_length, "%d", lengths[i]);
		
		char this_d_p_bbd[STRSIZE];
		sprintf(this_d_p_bbd, "%.2e", d_p_bbd[i]);
		
		char this_d_p_binomial[STRSIZE];
		sprintf(this_d_p_binomial, "%.2e", d_p_binomial[i]);
		
		char this_d_p_bbd_cor[STRSIZE];
		sprintf(this_d_p_bbd_cor, "%.2e", d_p_bbd_cor[i]);
		
		char this_d_p_binomial_cor[STRSIZE];
		sprintf(this_d_p_binomial_cor, "%.2e", d_p_binomial_cor[i]);
		
		char this_d_p_bbd_adj[STRSIZE];
		sprintf(this_d_p_bbd_adj, "%.2e", d_p_bbd_adj[i]);
		
		char this_d_p_bbd_cor_adj[STRSIZE];
		sprintf(this_d_p_bbd_cor_adj, "%.2e", d_p_bbd_cor_adj[i]);
		
		char this_d_p_binomial_adj[STRSIZE];
		sprintf(this_d_p_binomial_adj, "%.2e", d_p_binomial_adj[i]);
		
		char this_d_p_binomial_cor_adj[STRSIZE];
		sprintf(this_d_p_binomial_cor_adj, "%.2e", d_p_binomial_cor_adj[i]);
		
		// Round reptiming
		double rounded_reptiming = final_rounding(reptimings[i]);
		char this_reptiming[STRSIZE];
		sprintf(this_reptiming, "%.3f", rounded_reptiming);
		
		string outstring = annotation_mutation_counts[i][0] + "\t" +
											 annotation_mutation_counts[i][1] + "\t" +
											 annotation_mutation_counts[i][2] + "\t" +
											 annotation_mutation_counts[i][3] + "\t" +
											 annotation_mutation_counts[i][4] + "\t" +
											 annotation_mutation_counts[i][5] + "\t" +
											 string(this_reptiming) + "\t" +
											 string(this_length) + "\t" +
											 string(this_d_p_bbd) + "\t" +
											 string(this_d_p_binomial) + "\t" +
											 string(this_d_p_bbd_cor) + "\t" +
											 string(this_d_p_binomial_cor) + "\t" +
											 string(this_d_p_bbd_adj) + "\t" +
											 string(this_d_p_bbd_cor_adj) + "\t" +
											 string(this_d_p_binomial_adj) + "\t" +
											 string(this_d_p_binomial_cor_adj) + "\n";
		fprintf(outfile_ptr, outstring.c_str());
	}
	fclose(outfile_ptr);
	
	// Wrap up by removing the temporary files created along the way
	string rmcom = "rm " + mut_count;
	system(rmcom.c_str());
	rmcom = "rm " + gene_count;
	system(rmcom.c_str());
	rmcom = "rm " + pgene_count;
	system(rmcom.c_str());
	rmcom = "rm " + avg_infile;
	system(rmcom.c_str());
	rmcom = "rm " + avg_outfile;
	system(rmcom.c_str());
	
	// Reporting running time
	time_t t2;
	time(&t2);
	double run_sec = difftime(t2,t1);
	// Convert to minutes
	int run_min = (int)run_sec/60;
	if (run_min > 0) {
		run_sec = run_sec - ((double)run_min*60.0);
		printf("LARVA running time: %d min %d sec\n", run_min, (int)run_sec);
	} else {
		printf("LARVA running time: %d sec\n", (int)run_sec);
	}
	
	return 0;
}
