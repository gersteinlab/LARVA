# A custom-built CDF for the distribution of gene lengths.

expr_freqtable = read.table("expr_freqtable.txt")

pexpr <- function (x) {
	# Hardcode the frequency table creation
	# freqtable = read.table("expr_freqtable.txt")
	total = sum(expr_freqtable[,2])
	
	x = as.numeric(x)
	# print(expr_freqtable)
	# quit()
	
	# Area under the curve
	area = 0
	
	for (i in 1:nrow(expr_freqtable)) {
		# Sum up everything up to and including "x"
		
		# print(expr_freqtable[i,1])
		
		if (expr_freqtable[i,1] > x) {
			break
		} else {
			area = area + expr_freqtable[i,2]
		}
	}
	# print(area)
	# Do not want pexpr to be exactly equal to 1, or else (1-pexpr) will be 0
	if (area == total) {
		area = area-1
	}
	
	area/total
}

# arglist = commandArgs()
# print(arglist)
# preptime(arglist[6])
