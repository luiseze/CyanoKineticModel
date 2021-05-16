# concatenate the individual .tab files containing FCCs
import os

# directory with individual files
indir = 'parameter_sampling_results/fccs'

header_file = open('parameter_sampling_results/cbb_reaction_header.long.txt', 'r')

header = ['Conc_set', 'Stable_set', 'Reaction']

# concatenated fcc information (output file)
outfile = open('parameter_sampling_results/concset_stabstate_rxn_FCCs.tab', 'w')

for line in header_file.readlines():
	line = line.strip("\n")
	line = line.strip("\t")
	line = line.strip(" ")
	header.append(line)

outfile.write("\t".join(header) + "\n")

# concatenate individual fcc files in indir 
for filename in os.listdir(indir):
	file = open(os.path.join(indir, filename), 'r')
	for line in file.readlines():
		outfile.write(line)

outfile.close()