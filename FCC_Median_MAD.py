# read concset_stabstate_rxn_FCCs.tab and claculate median, MAD

import itertools
import pandas as pd
import numpy as np
from scipy import stats
import statistics

# read header line
fcc_infile = open('concset_stabstate_rxn_FCCs.tab', 'r')
header_line = fcc_infile.readline()
header_line = header_line.strip("\n")
header_list = header_line.split("\t")

fcc_infile.close()

# initialize output file
outfile = open('FCC_statistics_Median.tab', 'w')
outfile.write('Target' + "\t" + 'Effector' + "\t" + 'Median_FCC' + "\t" + "MAD" + "\n")

# 51 reactions, loop 51x through file to collect information for each reaction
# (reaction 1, reaction 2 ... reaction 51, reaction 1, reaction 2 ... reaction 51)
# Only read every 51st line in each loop
for i in range(1, 52):
	TargetNr = i
	
	# empty list for each effector reaction
	temp_dict = {}
	for rxns in range(3, len(header_list)):
		temp_dict[rxns] = []

	# store all FCCs in a list in a dictionary to later calculate median and median absolute deviation (MAD)
	with open('concset_stabstate_rxn_FCCs.tab', 'r') as fcc_infile:
		for line in itertools.islice(fcc_infile, i, None, 51):
			line = line.strip("\n")
			line_list = line.split("\t")
			rxn_counter = 0
			for val in line_list:
				if rxn_counter <3:
					rxn_counter += 1
				else:
					temp_dict[rxn_counter].append(float(val))
					rxn_counter += 1

	target = header_list[TargetNr+2]

	# calculate median and write it to output file
	for k, v in temp_dict.items():
		effector = header_list[k]
		FCCmedian = statistics.median(v)
		pd_v = pd.Series(v)
		FCCmad = stats.median_abs_deviation(pd_v)

		outfile.write(target + "\t" + effector + "\t" + str(FCCmedian) + "\t" + str(FCCmad) + "\n")

	# repeat for next reaction

outfile.close()

