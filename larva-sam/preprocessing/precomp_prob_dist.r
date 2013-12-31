# This R script takes the list of gene name and covariates, and for each gene,
# computes the special probability function used to determine how much "gravity"
# this gene has for random mutations. The values are scaled to sum to 1, then 
# scaled so that all terms are whole numbers (divide all by smallest term).
# These numbers are then used downstream in the Perl script for generating
# random variant distributions
# Rather than a return a value, the output values are printed to a file

# arglist[6] is the first argument
# arglist[7] is the second argument

source("powerlaw.R", chdir = TRUE)
source("reptime_cdf.r", chdir = TRUE)
source("length_cdf.r", chdir = TRUE)
source("expr_cdf.r", chdir = TRUE)

precompute_prob_dist <- function(input_file, output_file) {
	covar = read.table(input_file)
	
	# DEBUG
	# print(covar[1,1])
	# quit()
	
	expr = as.numeric(covar[,2])
	reptime = as.numeric(covar[,3])
	hic = as.numeric(covar[,4])
	length = as.numeric(covar[,5])
	
	# DEBUG
	# print(reptime)
	
	maxrow = nrow(covar)
	
	# DEBUG
	# print(maxrow)
	
	# The output matrix is "result"
	
	for (i in 1:maxrow) {
		this_prob = (1-pexpr(expr[i]))*
								preptime(reptime[i])*
								pnorm(hic[i], mean=20.0385567, sd=24.4418056, lower.tail=FALSE)*
								plength(length[i])
								# ppowerlaw(hic[i], alpha=1.49666501435019, xmin=4)
								
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
precompute_prob_dist(arglist[6], arglist[7])
