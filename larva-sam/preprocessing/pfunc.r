# A custom-built CDF for any of the distributions involved in the special CDF
# probability function
# Arguments:
# opt: Specify use of exome or whole genome distribution
# dist: Which distribution to use. Acceptable values: "reptime", "length", "expr"
# x: The quantile to evaluate for the given distribution

reptime_freqtable_ex = read.table("reptime_freqtable.txt")
reptime_freqtable_wg = read.table("reptime_freqtable_wg.txt")
lengths_freqtable = read.table("lengths_freqtable.txt")
expr_freqtable = read.table("expr_freqtable.txt")
h3k4me1_freqtable = read.table("h3k4me1_freqtable.txt")
h3k4me3_freqtable = read.table("h3k4me3_freqtable.txt")
expr_wg_freqtable = read.table("expr_wg_freqtable.txt")
snv_density_freqtable = read.table("snv_density_freqtable.txt")

pfunc <- function (opt, dist, x) {

	# So we ignore "opt" for all distributions except for "reptime"
	
	if (dist == "reptime" && opt == "e") { # Reptime exome version
		freqtable = reptime_freqtable_ex
	} else if (dist == "reptime" && opt == "g") { # Reptime whole genome version
		freqtable = reptime_freqtable_wg
	} else if (dist == "length") {
		freqtable = lengths_freqtable
	} else if (dist == "expr" && opt == "e") {
		freqtable = expr_freqtable
	} else if (dist == "expr" && opt == "g") {
		freqtable = expr_wg_freqtable
	} else if (dist == "H3K4me1") {
		freqtable = h3k4me1_freqtable
	} else if (dist == "H3K4me3") {
		freqtable = h3k4me3_freqtable
	} else if (dist == "snv_density") {
		freqtable = snv_density_freqtable
	} else { # This is not a valid combination of options
		print("Invalid input options")
		quit()
	}
	
	total = sum(freqtable[,2])
	
	x = as.numeric(x)
	# print(freqtable)
	# quit()
	
	# Area under the curve
	area = 0
	
	for (i in 1:nrow(freqtable)) {
		# Sum up everything up to and including "x"
		
		# print(freqtable[i,1])
		
		if (freqtable[i,1] > x) {
			break
		} else {
			area = area + freqtable[i,2]
		}
	}
	# print(area)
	
	# Do not want pfunc(dist="expr") (e and g) to be exactly equal to 1, or else (1-pfunc(dist="expr")) will be 0
	# Ditto the h3k4me1 and h3k4me3 distributions
	if ((dist == "expr" | dist == "H3K4me1" | dist == "H3K4me3") & area == total) {
		area = area-1
	}
	
	area/total
}
