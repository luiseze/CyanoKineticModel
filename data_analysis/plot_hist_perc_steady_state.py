# plot histogram for percent steady state

infile = open('parameter_sampling_results/met_set_vs_percent_steady.tab', 'r')

PercSteady = []

import matplotlib.pyplot as plt
import statistics
unstable_dict = {}
for line in infile.readlines():
	line = line.strip("\n")
	line_list = line.split("\t")
	PercSteady.append(float(line_list[1]))

# calculate min, max, median of the percentage of steady states
PercSteadyMin = min(PercSteady)
print('Min: ' + str(PercSteadyMin))
PercSteadyMax = max(PercSteady)
print('Max: ' + str(PercSteadyMax))
PercSteadyMed = statistics.median(PercSteady)
print('Median: ' + str(PercSteadyMed))

plt.rcParams.update({'font.size': 13})

plt.hist(PercSteady, color = "#8073ac", bins = 30, alpha = 0.4)
plt.xlabel('% stable steady states')
plt.ylabel('Number of sets')

plt.savefig("parameter_sampling_results/hist_perc_steady.pdf")