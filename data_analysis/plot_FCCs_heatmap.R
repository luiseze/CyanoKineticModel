# adapted from Janasch, M., Asplund-Samuelsson, J., Steuer, R., & Hudson, E. P. (2019). Kinetic modeling of the Calvin cycle identifies flux control and stable metabolomes in Synechocystis carbon fixation. Journal of Experimental Botany, 70(3), 1017–1031. https://doi.org/10.1093/jxb/ery382

fccs_file = "parameter_sampling_results/concset_stabstate_rxn_FCCs.tab.gz"
header_file = "parameter_sampling_results/cbb_reaction_header.long.txt"

# The enzyme in the column influences the flux of reactions in rows (FCC)

# Load custom labels for reactions
custom_rxn_labels = scan(header_file, character(), quote = "")

# Plot the influence distributions
library(ggplot2)
library(reshape2)

mean_mad_df = read.table("/Users/luise/Documents/Uni/Master/4_semester/Modeling/my_model/parameter_results_FBPase_5000x1000/FCC_statistics_Median.tab", header = TRUE)
fccs_med_mad = mean_mad_df

# Specify order of Effectors and Targets for plotting
fccs_med_mad$Effector = factor(fccs_med_mad$Effector, levels=custom_rxn_labels)
fccs_med_mad$Target = factor(fccs_med_mad$Target, levels=rev(custom_rxn_labels))

library(scales)

gp = ggplot(fccs_med_mad, aes(x=Target, y=Effector, alpha=MAD, fill=Median_FCC))
gp = gp + geom_tile(colour="white", size=0.5)
gp = gp + theme_bw()
gp = gp + scale_fill_gradientn(colours=c("#fdae61","#ffffbf","#2c7bb6"), values=rescale(c(min(fccs_med_mad$Median_FCC),0,max(fccs_med_mad$Median_FCC)), c(0,1)))
gp = gp + scale_alpha_continuous(range=c(1,0.1))
gp = gp + scale_x_discrete(position="top")
# gp = gp + theme(axis.text.x = element_text(angle = 90, hjust = 0, size = 7, face = "italic", vjust = 0.5))
gp = gp + theme(axis.text.x.top = element_text(vjust = 0.5, angle = 90, hjust = 0, size = 7, face = "italic"))
gp = gp + theme(axis.text.y = element_text(size = 7, face = "italic"))
gp = gp + theme(aspect.ratio=1)

outfile = paste(dirname(fccs_file), "FCCs_heatmap_Median.pdf", sep="/")
ggsave(outfile, gp, height=150/25.4, width=180/25.4)
