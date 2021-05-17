import sys
import random

# infile with sampled metabolite concentrations
infile = open(sys.argv[1], 'r')
outfile = open('Metabolite_Pool_Concentrations.txt', 'w')

header = True

for line in infile.readlines():
	line = line.strip('\n')
	
	# first line, identify indices for P_i, FUM, COA
	if header:
		# list of metabolites
		met_list = line.split('\t')
		header = False

		i = 0
		for met in met_list:
			if met == 'P_i':
				P_i_index = i
			elif met == 'FUM':
				FUM_index = i
			elif met == 'COA':
				COA_index = i
			else:
				pass
			i += 1

		outfile.write("\t".join(met_list) + "\t" + "PPool" + "\t" + "FUMPool" + "\t" + "COAPool" + "\n")

	# read metabolite concentrations, extract concentrations for P_i, FUM, COA and sample corresponding Pool concentration
	else:
		conc_list = line.split('\t')
		P_i = float(conc_list[P_i_index])
		FUM = float(conc_list[FUM_index])
		COA = float(conc_list[COA_index])

		P_i_Pool = P_i * random.uniform(1.1,5)
		FUM_Pool = FUM * random.uniform(1.1,5)
		COA_Pool = COA * random.uniform(1.1,5)

		outfile.write("\t".join(conc_list) + "\t" + str(P_i_Pool) + "\t" + str(FUM_Pool) + "\t" + str(COA_Pool) + "\n")
