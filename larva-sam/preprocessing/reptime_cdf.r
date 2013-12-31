# A custom-built CDF for the distribution of DNA replication times.

reptime_freqtable = read.table("reptime_freqtable.txt")

preptime <- function (x) {
	# Hardcode the frequency table creation
	# freqtable = read.table("reptime_freqtable.txt")
	total = sum(reptime_freqtable[,2])
	
	x = as.numeric(x)
	# print(reptime_freqtable)
	# quit()
	
	# Area under the curve
	area = 0
	
	for (i in 1:nrow(reptime_freqtable)) {
		# Sum up everything up to and including "x"
		
		# print(reptime_freqtable[i,1])
		
		if (reptime_freqtable[i,1] > x) {
			break
		} else {
			area = area + reptime_freqtable[i,2]
		}
	}
	# print(area)
	
	area/total
}

# arglist = commandArgs()
# print(arglist)
# preptime(arglist[6])
