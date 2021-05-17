# adapted from Janasch, M., Asplund-Samuelsson, J., Steuer, R., & Hudson, E. P. (2019). Kinetic modeling of the Calvin cycle identifies flux control and stable metabolomes in Synechocystis carbon fixation. Journal of Experimental Botany, 70(3), 1017–1031. https://doi.org/10.1093/jxb/ery382

conc_file = "Metabolite_Pool_Concentrations.txt"
perc_file = "met_set_vs_percent_steady.tab"
meta_file = "concentrations_and_cofactor_ratios_from_2016.tab"

library(reshape2)

# Load data
conc_data = read.table(conc_file, header=T, sep="\t")
perc_data = read.table(perc_file, header=F, sep="\t")
meta_data = read.table(meta_file, header=T, sep="\t")
colnames(perc_data) = c("set", "stable")

conc_data = conc_data[ , which(names(conc_data) %in% c("FBP", "RuBP"))]

meta_data = meta_data[grep("RuBP", meta_data$metabolite),]

# Transform the data
conc_data = log10(conc_data)
meta_data$concentration = log10(meta_data$concentration)

# Use Synechocystis and Synechococcus data
meta_data = subset(meta_data, Organism %in% c("PCC 6803", "PCC 7942"))

# Split data into lower and upper stable steady state categories
perc_data = perc_data[order(perc_data$stable, decreasing=T),]

# Divide metabolite concentration sets into high and low % stable steady states

# By Deciles
deciles = quantile(perc_data$stable, prob = seq(0, 1, length = 11), type = 5)
d10 = deciles[10]
d1 = deciles[2]

perc_hi = subset(perc_data, stable > d10)
perc_lo = subset(perc_data, stable <= d1)

# Subset concentrations to those in high and low % stable steady states groups
conc_hi = conc_data[perc_hi$set,]
conc_lo = conc_data[perc_lo$set,]

# Plot it
library(ggplot2)

c_lo_long = reshape2::melt(conc_lo)
colnames(c_lo_long) = c("metabolite", "concentration")
c_lo_long$stability = rep("unstable", nrow(c_lo_long))

c_hi_long = reshape2::melt(conc_hi)
colnames(c_hi_long) = c("metabolite", "concentration")
c_hi_long$stability = rep("stable", nrow(c_hi_long))

c_long = rbind(c_hi_long, c_lo_long)

meta_data$stability = "stable"

# Plot metabolite concentrations separately

gp = ggplot(
  c_long[grep("(charge)|(ratio)", c_long$metabolite, perl=T, invert=T),],
  aes(x=concentration, fill=stability)
)
gp = gp + geom_density(alpha=0.4)
gp = gp + geom_point(
  aes(x=concentration, colour=Organism, shape=Organism),
  subset(meta_data, metabolite %in% grep("(charge)|(ratio)", c_long$metabolite, perl=T, invert=T, value=T)),
  y=0, show.legend=F, colour="white", alpha=0.7, size=2, shape=19
)
gp = gp + geom_point(
  aes(x=concentration, colour=Organism, shape=Organism),
  subset(meta_data, metabolite %in% grep("(charge)|(ratio)", c_long$metabolite, perl=T, invert=T, value=T)),
  y=0, show.legend=F, fill="black", alpha=0.8, size=1, stroke=0.6
)
gp = gp + scale_shape_manual(values=c(6, 2))
gp = gp + scale_colour_manual(values=c("black", "red"))
gp = gp + theme_bw()
gp = gp + facet_wrap(~metabolite, ncol=6)
gp = gp + scale_fill_manual(values=c("#8073ac","#e08214"))
gp = gp + theme(strip.background = element_blank())

gp = gp + xlim(-4,2)
gp = gp + ylim(0,2.5)

outfile = paste(dirname(perc_file), "concs_vs_stability.with_meta_testtest.pdf", sep="/")
ggsave(outfile, gp, width=200/25.4, height=50/25.4)
