# Import libraries
import numpy as np
import pandas as pd

########################

# from SMatrix_from_reactions import read_equations
# from SMatrix_from_reactions import create_Smatrix

# reactions_infile = open('reactions.txt', 'r')

# reaction_dict, seen_met = read_equations(reactions_infile)
# S_pandas, S = create_Smatrix(reaction_dict, seen_met)

########################

conc_infile = open('concentration_ranges.txt', 'r')

def read_ranges(conc_infile):
	# read allowed concentration ranges and ratios between different metabolites from the input file

	conc_dict = {}
	ratio_dict = {}
	met_ratio_dict = {} # with keys being metabolites and values being the corresponding ratios (as seen in ratio_dict.keys())

	for lines in conc_infile.readlines():

		# header line and pool ranges
		if lines.startswith('NAME') or 'Pool' in lines:
			pass

		# ratios between different metabolites
		elif '/' in lines:
			lines = lines.strip('\n')
			dividend = lines.split('\t')[0].split('/')[0]
			divisor = lines.split('\t')[0].split('/')[1]
			ratio_dict[lines.split('\t')[0]] = {'dividend': dividend, 'divisor': divisor, 'min': lines.split('\t')[1], 'max': lines.split('\t')[2]}
			met_ratio_dict[dividend] = lines.split('\t')[0]
			met_ratio_dict[divisor] = lines.split('\t')[0]

		# ranges for metabolites
		else:
			lines = lines.strip('\n')
			conc_dict[lines.split('\t')[0]] = {'min': lines.split('\t')[1], 'max': lines.split('\t')[2]}

	return conc_dict, ratio_dict, met_ratio_dict

def ranges_to_array(conc_dict, S_pandas):
	# convert the ranges to numpy array

	for i in range(0, len(S_pandas.index)):
		met = S_pandas.index[i]
		minimum = float(conc_dict[met]['min'])/1000
		maximum = float(conc_dict[met]['max'])/1000

		if 'conc_array' in locals():
			new_conc = np.array([minimum, maximum])
			conc_array = np.vstack((conc_array, new_conc))

		else:
			conc_array = np.array([minimum, maximum])

	return conc_array

def ratios_to_array(ratio_dict, met_ratio_dict, S_pandas):
	# convert the ratios to numpy array

	seen_ratios = []

	S_position = {}
	for i in range(0, len(S_pandas.index)):
		met = S_pandas.index[i]
		S_position[met] = int(i)

	for i in range(0, len(S_pandas.index)):
		met = S_pandas.index[i]
		# create numpy array including the given ratios
		if met in met_ratio_dict.keys():
			if met_ratio_dict[met] not in seen_ratios:
				seen_ratios.append(met_ratio_dict[met]) 
				minimum = float(ratio_dict[met_ratio_dict[met]]['min'])
				maximum = float(ratio_dict[met_ratio_dict[met]]['max'])
				
				if 'ratio_limits' in locals():
					new_ratio = np.array([minimum, maximum])
					ratio_limits = np.vstack((ratio_limits, new_ratio))
				else:
					ratio_limits = np.array([minimum, maximum])
		
				# create numpy array including for which metabolites the ratios apply

				dividend = ratio_dict[met_ratio_dict[met]]['dividend']
				divisor = ratio_dict[met_ratio_dict[met]]['divisor']

				pos_dividend = S_position[dividend]
				pos_divisor = S_position[divisor]

				new_vector = np.zeros((len(S_pandas.index), 1))
				new_vector[pos_dividend, 0] = 1
				new_vector[pos_divisor, 0] = -1

				if 'ratio_mat' in locals():
					ratio_mat = np.hstack((ratio_mat, new_vector))
				else:
					ratio_mat = new_vector

			else:
				pass
		else:
			pass

	return ratio_limits, ratio_mat
