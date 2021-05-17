# read delta G values from equilibrator_results.tsv

infile = open('equilibrator_results.tsv', 'r')

import numpy as np

def read_dg(infile):
	dg_list = []
	for line in infile:
		if not line.startswith("'"): # skips line with headers
			pass
		else:
			line = line.strip("\n")
			line_list = line.split("\t")
			dg = line_list[2].split(" ")[0]
			dg = dg.strip("(")
			dg_list.append(float(dg))

	dg_array = np.array(dg_list)

	return(dg_array)
