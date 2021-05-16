# Import libraries
import numpy as np
import re
import pandas as pd

# opening infile: .txt file created from excel sheet containing the reactions in the following format: 'A + B => C' or 'A + B <=> C'
infile = open('reactions.txt', 'r')

# function reading the equations from the infile and storing information in dictionaries
def read_equations(infile):

	reaction_id = 1 # index of reaction, starting from 1

	reaction_dict = {} # nested dict, key = reaction_id, value = dict with key = 'reactants'/'products' and value = list of reactants/products 

	seen_met = {} # keep track of metabolites, key = metabolite, value = list of reaction_ids

	for line in infile: # reads reactions from file (one reaction per line in format 'A + B => C' or 'A + B <=> C')
		line = line.strip('\n')
		line = line.strip("'")
		met_list = line.split("+") # splits reactions based on '+'
		side = 'reactants' # starting with the left side of the equation
		# initializing dictionary as described above
		reaction_dict[reaction_id] = {}
		reaction_dict[reaction_id]['reactants'] = []
		reaction_dict[reaction_id]['products'] = []
		# save information from met_list in dictionary
		for met in met_list:
			met = met.strip(' ')
			if '=>' not in met and '<=>' not in met: # if it is a metabolite, add it to dictionary
				if ' ' in met: # if stoichiometry is given for the metabolite such as 2 P3G, then it is split into '2' and 'P3G'
					reaction_dict[reaction_id][side].append((met.split(' ')[1], int(met.split(' ')[0])))
					met = met.split(' ')[1]
				else: # if no stoichiometry is given, the tuple is (metabolite, 1) - 1 for 1x metabolite
					reaction_dict[reaction_id][side].append((met, 1))
				# save the metabolites in seen_met list with format as described above
				if met not in seen_met.keys(): 
					seen_met[met] = [reaction_id]
				else:
					seen_met[met].append(reaction_id)
			elif '=>' in met or '<=>' in met: # in this case switch from reactants to products (note: met in this case includes the last metabolite before and first metabolite after => or <=>)
				met_list2 = re.split('<=>|=>', met) # entry which contains last reactant and first product; split it into reactant and product
				for i in range(0, 2): # the first one will be added as reactant, the second one as product
					met = met_list2[i].strip(' ')

					if i == 0:
						side = 'reactants'
					elif i == 1:
						side = 'products'

					if ' ' in met: # if stoichiometry is given
						reaction_dict[reaction_id][side].append((met.split(' ')[1], met.split(' ')[0]))
						met = met.split(' ')[1]
					else:
						reaction_dict[reaction_id][side].append((met, 1))

					if met not in seen_met.keys(): 
						seen_met[met] = [reaction_id]
					else:
						seen_met[met].append(reaction_id)

		reaction_id += 1

	return reaction_dict, seen_met

def create_Smatrix(reaction_dict, seen_met):
	
	num_reactions = len(reaction_dict) # the number of equations equals the number of keys in the reaction_dict dictionary
	num_metabolites = len(seen_met) # the number of metabolites equals the number of keys in the seen_met dictionary

	S = pd.DataFrame(np.zeros((num_metabolites, num_reactions))) # create S matrix filled with zeros for now

	# each metabolite is assigned one metabolite id starting with 0, information stored in this dictionary: keys = metabolites, values = metabolite ids
	met_id = {}

	# metabolite id is assigned as described above
	metabolite_id = 0
	for met in seen_met.keys():
		met_id[met] = metabolite_id
		metabolite_id += 1

	S.columns = reaction_dict.keys() # column names now from 1 to ... (number of reactions)
	S.index = seen_met.keys() # row names of the data frame are the names of the metabolites

	for k,v in reaction_dict.items(): # k is reaction_id, v is another dictionary which is looped through in the following for loop
		for key, value in v.items(): # keys are 'reactants' and 'products' and values are a list for each key containing the reactants and products of the reaction, respectively.
			# -1 for reactants and 1 for products in S matrix, in following for loop multiplied with stoichiometry
			if key == 'reactants':
				entry = -1
			if key == 'products':
				entry = 1

			# filling S matrix with information abour reactions, S[k] is the column of the reaction, S[k][m[0]] selects the row of a certain metabolite (m[0])
			for m in value: # loop through metabolites (tuples with (metabolite, stoichiometry), m[1] is stoichiometry)
				S[k][m[0]] = float(entry * m[1])
				pass

	# convert pandas data frame to numpy array as this is the format needed in the hin and run metabolite sampling
	S_numpy = S.to_numpy()

	return S, S_numpy
