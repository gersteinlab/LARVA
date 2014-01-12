# This R script takes the list of gene name and covariates, and for each gene,
# computes the special probability function used to determine how much "gravity"
# this gene has for random mutations. The values are scaled to sum to 1, then 
# scaled so that all terms are whole numbers (divide all by smallest term).
# These numbers are then used downstream in the Perl script for generating
# random variant distributions
# Rather than a return a value, the output values are printed to a file

# arglist[6] is the first argument
# arglist[7] is the second argument
# arglist[8] is the third argument

# Arguments:
# opt: Either "e" for exome, or "g" for whole genome
# input_file: The input file contains the names of all the regions in the first
# column, and an additional column for each of the attributes' values
# output_file: The output file contains the names of all the regions in the first
# column, and a second column with the value assigned by the CDF probability function

source("pfunc.r", chdir = TRUE)

precompute_prob_dist <- function(opt, input_file, output_file) {
	covar = read.table(input_file)
	
	# DEBUG
	# print(covar[1,1])
	# quit()
	
	if (opt == "e") { # Exome version
		expr = as.numeric(covar[,2])
		reptime = as.numeric(covar[,3])
		hic = as.numeric(covar[,4])
		length = as.numeric(covar[,5])
	} else if (opt == "g") { # Whole genome version
		reptime = as.numeric(covar[,2])
		h3k4me1 = as.numeric(covar[,3])
		h3k4me3 = as.numeric(covar[,4])
		expr = as.numeric(covar[,5])
		snv_density = as.numeric(covar[,6])
	} else { # This is not a valid option
		print("Option must be either \'e\' for generating exome probability distribution, or \'g\' for generating whole genome probability distribution")
		quit()
	}
	
	# DEBUG
	# print(reptime)
	
	maxrow = nrow(covar)
	
	# DEBUG
	# print(maxrow)
	
	# The output matrix is "result"
	
	for (i in 1:maxrow) {
		if (opt == "e") { # Exome version
			this_prob = log(1-pfunc("e", "expr", expr[i]))+
									log(pfunc("e", "reptime", reptime[i]))+
									log(1-pnorm(hic[i], mean=20.0385567, sd=24.4418056, lower.tail=FALSE))+
									log(pfunc("e", "length", length[i]))
		} else if (opt == "g") { # Whole genome version
			this_prob = log(pfunc("g", "reptime", reptime[i]))+
									log(1-pfunc("g", "H3K4me1", h3k4me1[i]))+
									log(1-pfunc("g", "H3K4me3", h3k4me3[i]))+
									log(1-pfunc("g", "expr", expr[i]))+
									log(pfunc("g", "snv_density", snv_density[i]))
		}
								
		# DEBUG
		# print(i)
		# print(pexpr(expr[i]))
		# print(preptime(reptime[i]))
		# print(pnorm(hic[i], mean=20.0385567, sd=24.4418056, lower.tail=FALSE))
		# print(plength(length[i]))
		# print(this_prob)
								
		if (exists("result")) {
			result = rbind(result, c(toString(covar[i,1]),this_prob))
		} else {
			result = rbind(c(toString(covar[i,1]),this_prob))
		}
		
		# DEBUG
		# this_sum = sum(as.numeric(result[,2]))
		# print(this_sum)
	}
	
	# result[,2] = as.numeric(result[,2])
	
	# DEBUG
	# print(result)
	
	# All probabilities computed, time for scaling
	# First, divide all values by the current sum
	first_sum = sum(as.numeric(result[,2]))
	result[,2] = as.numeric(result[,2])/first_sum
	
	# Now, divide all terms by the smallest term
	result[,2] = as.numeric(result[,2])/min(as.numeric(result[,2]))
	
	# Round the terms
	result[,2] = round(as.numeric(result[,2]))
	
	# OK, let us print this
	write.table(result, file=output_file, sep="\t")
}

arglist = commandArgs()
# print(arglist)
precompute_prob_dist(arglist[6], arglist[7], arglist[8])
