import math

from equilibrator_api import ComponentContribution, Q_
cc = ComponentContribution()
cc.p_h = Q_(8.4) # set pH = 8.4
cc.p_mg = Q_(3.0) # set pMg = 3.0
cc.ionic_strength = Q_("0.1M") # set ionic strength = 0.1M

from SMatrix_from_reactions import read_equations

reactions_infile = open('reactions.txt', 'r')
id_infile = open('kegg_ids.txt', 'r')
outfile = open('equilibrator_results.tsv', 'w')

reaction_dict, seen_met = read_equations(reactions_infile)

def read_ids(id_infile):
	id_dict = {} # dict with model_id as keys, values are further dictionaries with the type of id as keys and respective ids as values

	id_types = {}
	
	for line in id_infile:
		line = line.strip("\n")
		ids_list = line.split("\t")

		if line.startswith('model_id'):
			counter = 0
			for ids in ids_list:
				id_types[counter] = ids
				counter += 1

		else:
			id_dict[ids_list[0]] = {}
			for i in range(1, len(id_types.keys())):
				id_dict[ids_list[0]][id_types[i]] = ids_list[i]


	return(id_dict)

id_dict = read_ids(id_infile)

outfile.write("Reaction" + "\t" + "K_eq" + "\t" + "ΔG'0" + "\n")

counter = 1

R = 8.314 # [J/(mol*K)] = [kg*m2/(s2*mol*K)]
T = 303.15 # [K]

RT = R*T/1000 # /1000 because ΔG'0 is in kJ/mol

# needs to be opened again
reactions_infile = open('reactions.txt', 'r')

# convert reaction equations into suitable format using kegg ids
for line in reactions_infile:
	line = line.strip("\n")
	outfile.write(line)

	formula = str()

	for reactants in reaction_dict[counter]['reactants']:
		if reactants[1] == 1:
			formula += 'kegg:' + id_dict[reactants[0]]['kegg_id'] + " + "
		else:
			formula += reactants[1] + ' kegg:' + id_dict[reactants[0]]['kegg_id'] + " + "

	formula = formula.strip(" + ")
	formula += " = "

	for products in reaction_dict[counter]['products']:
		if products[1] == 1:
			formula += 'kegg:' + id_dict[products[0]]['kegg_id'] + " + "
		else:
			formula += products[1] + ' kegg:' + id_dict[products[0]]['kegg_id'] + " + "

	formula = formula.strip(" + ")

	# calculate ΔG'0 in kJ/mol, K_eq and write to output file
	rxn = cc.parse_reaction_formula(formula)
	std_dg_prime = cc.standard_dg_prime(rxn) # ΔG'0 in kJ/mol
	G0 = str(std_dg_prime).strip("(").split(" ")[0]
	K_eq = math.exp(-float(G0)/RT)
	K_eq = round(K_eq,5)

	outfile.write("\t" + str(K_eq) + "\t" + str(std_dg_prime) + "\n")

	counter += 1
