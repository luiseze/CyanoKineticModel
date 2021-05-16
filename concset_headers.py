# add headers to concset_stability_parameters_noheaders.tab

infile = open('parameter_sampling_results/concset_stability_parameters_noheaders.tab','r')

header_file = open('parameter_sampling_results/par/par_header.long.txt', 'r')

outfile = open('parameter_sampling_results/concset_stability_parameters.tab','w')

header = ['Conc_set', 'Stable']

for line in header_file.readlines():
	line = line.strip("\n")
	line = line.strip("\t")
	line = line.strip(" ")
	header.append(line)

outfile.write("\t".join(header) + "\n")

for line in infile.readlines():
	outfile.write(line)

outfile.close()