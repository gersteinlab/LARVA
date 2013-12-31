# A custom-built CDF for the distribution of gene lengths.

lengths_freqtable = read.table("lengths_freqtable.txt")

plength <- function (x) {
	# Hardcode the frequency table creation
	# freqtable = read.table("lengths_freqtable.txt")
	total = sum(lengths_freqtable[,2])
	
	x = as.numeric(x)
	# print(lengths_freqtable)
	# quit()
	
	# Area under the curve
	area = 0
	
	for (i in 1:nrow(lengths_freqtable)) {
		# Sum up everything up to and including "x"
		
		# print(lengths_freqtable[i,1])
		
		if (lengths_freqtable[i,1] > x) {
			break
		} else {
			area = area + lengths_freqtable[i,2]
		}
	}
	# print(area)
	
	area/total
}

# arglist = commandArgs()
# print(arglist)
# preptime(arglist[6])
