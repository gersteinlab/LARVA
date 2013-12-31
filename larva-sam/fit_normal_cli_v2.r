# Script to take a distribution, fit a Normal, and produce the 95% CI
# This version is intended to take its arguments from the command line
# Arguments:
# arglist[1]/input_file: The file with the vector of numbers to use as input
# arglist[2]/actual_data: This number represents the observed data, for p-value determination
# arglist[3]/output_file: The file to which the function's output is to be written
# v2 adjusts the script to return a result vector, rather than write to a file

# global_result/local_result layout:
# [1] Mean
# [2] SD
# [3] 95% CI lower bound
# [4] 95% CI upper bound
# [5] p-value

library(MASS)

fit_normal <- function(input_file, actual_data) {

	invec = read.table(input_file)
	invec = as.numeric(invec)
	
	fit = fitdistr(invec, "normal")
	# write.table(fit, file=output_file, sep = "\t")
	
	# DEBUG
	# print("Begin write")
	
	# write("Distribution stats: ", file=output_file, append = FALSE, sep = "\t")
	# write("Mean: ", file=output_file, append = TRUE, sep = "\t")
	local_result = c(fit$estimate[1])
	# write("SD: ", file=output_file, append = TRUE, sep = "\t")
	local_result = c(local_result, fit$estimate[2])
	
	# write("95% CI: ", file=output_file, append = TRUE, sep = "\t")
	local_result = c(local_result, fit$estimate[1]-(2*fit$estimate[2]))
	local_result = c(local_result, fit$estimate[1]+(2*fit$estimate[2]))
	
	pvalue = 0
	
	# write("Actual data: ", file=output_file, append = TRUE, sep = "\t")
	# write(actual_data, file=output_file, append = TRUE, sep = "\t")
	
	if (as.numeric(actual_data) < fit$estimate[1]) {
		# write("Depletion", file=output_file, append = TRUE, sep = "\t")
		pvalue = pnorm(as.numeric(actual_data), mean = fit$estimate[1], sd = fit$estimate[2], lower.tail = TRUE)
	} else {
		# write("Enrichment", file=output_file, append = TRUE, sep = "\t")
		pvalue = pnorm(as.numeric(actual_data), mean = fit$estimate[1], sd = fit$estimate[2], lower.tail = FALSE)
	}
	
	# write("p-value: ", file=output_file, append = TRUE, sep = "\t")
	local_result = c(local_result, pvalue)
	return(local_result)
}

arglist = commandArgs()
# print(arglist)
global_result = fit_normal(arglist[6], arglist[7])
global_result
