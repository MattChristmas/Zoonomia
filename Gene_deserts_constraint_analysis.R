library(readr)
library(ggplot2)
library(boot)
library(ggrepel)
library(dplyr)
library(broom)
library(vegan)
library(SoDA)
library(ggpubr)
library(gridExtra)
library(cowplot)


# This script is for analysing constraint in gene deserts near developmental transcription factors
# Developmental transcription factors are those used by the following publication:
#Touceda-Suárez, María, et al. "Ancient Genomic Regulatory Blocks Are a Source for Regulatory 
#Gene Deserts in Vertebrates after Whole-Genome Duplications." Molecular biology and evolution 
#37.10 (2020): 2857-2864.

options(scipen=10000)
# A function to calculate mean from the data for bootstrapping

samplemean <- function(x, d) {  # This is the function to calculate mean from the data (x), using a bootstrap sample, d.
  return(mean(x[d]))
}

# function for number of observations
give.n <- function(x){
  return(c(y = median(x)*1.2, label = length(x))) 
  # experiment with the multiplier to find the perfect position
}

# Function for calculating the standard error
se <- function(x) sqrt(var(x)/length(x))

# import intergenic phylop file

intergen <- read_delim("ensembl_intergenic_regions_phylop_td_genes_defined.txt", 
                          "\t", escape_double = FALSE, trim_ws = TRUE)
# Remove the few massive regions (>4mil) - all are non-td
intergen <- subset(intergen, intergenic_distance < 4000000)

# And those below 0
intergen <- subset(intergen, intergenic_distance > 0)

# Import list of td-genes that are in ancient genomic regulatory blocks (GRBs)

grb_genes <- scan("GRB_genes.txt",what = list(name = character()))

# And use this list to identify genes in main df as GRB or not
intergen$GRB = "FALSE"
find_grb_genes <- intergen$td_gene_name %in% grb_genes$name
intergen$GRB[find_grb_genes] = "TRUE"
table(intergen$GRB)

# plot boxplot of constraint for gene categories

intergen_boxplot <- ggplot(data=intergen, aes(x=GRB, y=prop_constraint)) +
  geom_boxplot() +
  xlab("Intergenic region") + 
  ylab("Proportion of constraint") +
  #scale_colour_manual(values=c("red","grey","orange")) +
  #scale_y_continuous(limits = c(0, 4.8), expand = c(0, 0)) +
  theme_bw() + 
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16,face="bold")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))

intergen_boxplot

# Calculate constraint means
td_bs <- intergen$category == "td_bs"
td_neighbour <- intergen$category == "td_neighbour"
non_td <- intergen$category == "non_td"
GRB <- intergen$GRB =="TRUE"

td_bs_mean_constraint <- mean(intergen$prop_constraint[td_bs])
td_neighbour_mean_constraint <- mean(intergen$prop_constraint[td_neighbour])
non_td_mean_constraint <- mean(intergen$prop_constraint[non_td])
GRB_mean_constraint <- mean(intergen$prop_constraint[GRB])

# Put td_bs and td_neighbour regions into same category
intergen$td_region = "No"
intergen$td_region[td_bs] = "Yes"
intergen$td_region[td_neighbour] = "Yes"
intergen$td_region[GRB] = "GRB"

# Look at size distribution of intergenic regions
td_only <- subset(intergen, category == "td_bs" | category == "td_neighbour")

intergen_size_distribution <- ggplot(data=td_only, aes(x=intergenic_distance/1000, fill=category)) +
  geom_histogram(binwidth=1,position="identity", alpha=0.5) +
  xlab("Intergenic distance (Kbp)") + 
  ylab("Frequency") +
  #scale_colour_manual(values=c("red","grey","orange")) +
  scale_x_continuous(limits = c(0, 500)) +
  theme_bw() + 
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16,face="bold")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))

intergen_size_distribution

#Look at the top 10% of intergenic length
quantile(intergen$intergenic_distance, 0.9)
# 200285.8 
top_10pc_intergenic_length <- subset(intergen, intergenic_distance>200000)

# Number of genes in the top 10% of intergenic lengths:
length(with(top_10pc_intergenic_length, union(gene_1_ID, gene_2_ID))) # 2,882 genes

# and how many of these are td-genes?
top_10pc_td_gene_subset <- subset(top_10pc_intergenic_length, category == "td_bs" | category == "td_neighbour")
length(unique(top_10pc_td_gene_subset$td_gene_ID)) # 374 genes

# And how many are in GRBs
top_10pc_GRB_subset <- subset(top_10pc_intergenic_length, GRB == "TRUE")
length(unique(top_10pc_GRB_subset$td_gene_ID)) # 51 genes

# Plot length of intergenic region against constraint

intergen_length_v_constraint <- ggplot(data=top_10pc_intergenic_length, aes(x=intergenic_distance/1000, y=prop_constraint, colour=td_region)) +
  geom_point(size=0.5) +
  #geom_point(size=0.8, colour = "blue", data=top_10pc_GRB_subset) +
  #geom_text_repel(aes(label=ifelse((prop_constraint>0.07 & (intergenic_distance/1000) > 1000), gene_1,""))) +
  geom_text_repel(aes(label=ifelse(prop_constraint>0.07 & (category == "td_bs" | category == "td_neighbour"), td_gene_name,"")), show.legend = FALSE) +
  xlab("Intergenic length (Kbp)") + 
  ylab("Proportion of constraint") +
  geom_vline(xintercept=200, lty=3, colour="tomato") + #95th pc of length
  geom_vline(xintercept=400, lty=3, colour="purple") + #95th pc of length
  scale_colour_manual(values=c("tomato","orange","grey"), name="Category", breaks = c("GRB","Yes", "No"),labels=c("GRB","trans-dev gene", "other")) +
  scale_x_continuous(breaks=c(200,400,1000,2000,3000,4000)) +
  theme_bw() + 
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16,face="bold")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))

intergen_length_v_constraint

# plot boxplot of constraint for gene categories

top_10pc_intergenic_length$td_region <- factor(top_10pc_intergenic_length$td_region,levels = c("GRB","Yes","No"),ordered = TRUE)

intergen_top10_boxplot <- ggplot(data=top_10pc_intergenic_length, aes(x=td_region, y=prop_constraint, fill=td_region)) +
  geom_boxplot() +
  xlab("Category") + 
  ylab("Proportion of constraint") +
  stat_summary(fun.data = give.n, geom = "text", fun = median) +
  scale_fill_manual(values=c("tomato","orange","grey"), name="Category", breaks = c("GRB","Yes", "No"),labels=c("GRB","trans-dev gene", "other")) +
  scale_x_discrete(breaks=c("GRB","Yes","No"), labels=c("GRB","trans-dev gene", "other")) +
  #scale_y_continuous(limits = c(0, 4.8), expand = c(0, 0)) +
  theme_bw() + 
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16,face="bold")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), legend.position = "none")

intergen_top10_boxplot

#and the top 5% of intergenic length
quantile(intergen$intergenic_distance, 0.95)
top_5pc_intergenic_length <- subset(intergen, intergenic_distance>400000)

# Remove GRB label
GRBs <- top_5pc_intergenic_length$GRB =="TRUE"
top_5pc_intergenic_length$td_region[GRBs] <- "Yes"

# plot boxplot of constraint for gene categories

top_5pc_intergenic_length$td_region <- factor(top_5pc_intergenic_length$td_region,levels = c("Yes","No"),ordered = TRUE)

intergen_top5_boxplot <- ggplot(data=top_5pc_intergenic_length, aes(x=td_region, y=prop_constraint, fill=td_region)) +
  geom_boxplot() +
  xlab("Developmental transcription factor neighbour?") + 
  ylab("Proportion of constraint") +
  stat_summary(fun.data = give.n, geom = "text", fun = median) +
  scale_fill_manual(values=c("tomato","grey"), name="Developmental transcription factor neighbour?", breaks = c("Yes", "No"),labels=c("Yes", "No")) +
  scale_x_discrete(breaks=c("Yes","No"), labels=c("Yes", "No")) +
  #scale_y_continuous(limits = c(0, 4.8), expand = c(0, 0)) +
  theme_bw() + 
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16,face="bold")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), legend.position = "none")

intergen_top5_boxplot

# Calculate constraint means for top 5 percentile of intergenic regions
#top5_GRB <- top_5pc_intergenic_length$td_region == "GRB"
top5_td <- top_5pc_intergenic_length$td_region == "Yes"
top5_non_td <- top_5pc_intergenic_length$td_region == "No"

#top5_GRB_mean_constraint <- mean(top_5pc_intergenic_length$prop_constraint[top5_GRB])
top5_td_mean_constraint <- mean(top_5pc_intergenic_length$prop_constraint[top5_td])
top5_non_td_mean_constraint <- mean(top_5pc_intergenic_length$prop_constraint[top5_non_td])

constraint_means <- data.frame(row.names = c(1,2))
constraint_means$type = c("trans-dev gene", "other")
constraint_means$mean = c(top5_td_mean_constraint, top5_non_td_mean_constraint)



# Calculate 95% confidence intervals using bootstrapping

boot_out = list() # Create lists for storing bootstrap outputs
boot_ci_out = list()

cats=c("Yes","No")

# Loop over each category and calculate bootstrap means and CIs of constraint

for (cat in 1:length(cats)) {
  
  temp=top_5pc_intergenic_length[top_5pc_intergenic_length$td_region==cats[cat],]
  boot_out[[cat]] <- boot(temp$prop_constraint, samplemean, R=1000)
  boot_ci_out[[cat]] <- boot.ci(boot_out[[cat]], type = "perc")
  # Extract confidence intervals and add to kin_means df
  constraint_means$lower.ci[cat] <- (boot_ci_out[[cat]][["percent"]][(4)])
  constraint_means$upper.ci[cat] <- (boot_ci_out[[cat]][["percent"]][(5)])
  
  #  i=i+1
}

# Set the order of categories for plotting
constraint_means$type = factor(constraint_means$type, levels = c("trans-dev gene", "other"))
# And plot
intergen_top5_mean_constraint <- ggplot(data=constraint_means, aes(x = type, y = mean)) + 
  geom_col(width = 0.5, fill="lightgrey") +
  geom_errorbar(ymin=constraint_means$lower.ci, ymax=constraint_means$upper.ci, width=0.2) +
  scale_y_continuous(expand=c(0,0),limits=c(0,0.04)) +
  xlab("") +
  ylab("Proportion of constraint") +
  scale_x_discrete(breaks=c("Yes","No"), labels=c("Yes", "No")) +
  scale_color_manual(values=c("grey","orange")) +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line=element_line(colour = "black"))
intergen_top5_mean_constraint

# Test significance of difference
t.test(top_5pc_intergenic_length$prop_constraint[top5_td],top_5pc_intergenic_length$prop_constraint[top5_non_td])
wilcox.test(top_5pc_intergenic_length$prop_constraint[top5_td],top_5pc_intergenic_length$prop_constraint[top5_non_td])
length(top_5pc_intergenic_length$prop_constraint[top5_td])
length(top_5pc_intergenic_length$prop_constraint[top5_non_td])
# Is there enrichment for trans-dev genes in the top 5% intergenic regions?

top_5pc_td_intergen = sum(top_5pc_intergenic_length$td_region == "GRB") + sum(top_5pc_intergenic_length$td_region == "Yes") 
top_5pc_non_intergen = sum(top_5pc_intergenic_length$td_region == "No")
td_intergen = sum(intergen$td_region == "GRB") + sum(intergen$td_region == "Yes") #3243
non_td_intergen = sum(intergen$td_region == "No")

# Create 2x2 matrix
intergen_type_contingency_table <- rbind(c(non_td_intergen,td_intergen),c(top_5pc_non_intergen,top_5pc_td_intergen))
rownames(intergen_type_contingency_table) <- c("total","top_5pc")
colnames(intergen_type_contingency_table) <- c("non_td","td")
fisher.test(intergen_type_contingency_table)
chisq.test(intergen_type_contingency_table)


## Look at enhancer constraint in intergenic regions

enhancers <- read_delim("encode_human_p_and_dELS_phyloP_TD_associations.txt", 
                       "\t", escape_double = FALSE, trim_ws = TRUE)

td_enhancers <- enhancers$intergenic_category == "td_neighbour" | enhancers$intergenic_category == "td_bs"
enhancers$TD_associated = "FALSE"
enhancers$TD_associated[td_enhancers] = "TRUE"

enhancer_constraint_plot <- ggplot(data=enhancers, aes(x = TD_associated, y = prop_constraint)) + 
  geom_boxplot() +
  #scale_y_continuous(limits=c(0,0.055)) +
  #xlab("Intergenic category") +
  #ylab("Proportion of constraint") +
  #scale_color_manual(values=c("grey","orange")) +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line=element_line(colour = "black"))
enhancer_constraint_plot

# Subset for the longest intergenic regions
enhancers_large_intergenic_regions <- subset(enhancers, intergenic_length >400000)

enhancer_constraint_large_intergenic_regions_plot <- ggplot(data=enhancers_large_intergenic_regions, aes(x = TD_associated, y = prop_constraint)) + 
  geom_boxplot() +
  #scale_y_continuous(limits=c(0,0.055)) +
  #xlab("Intergenic category") +
  #ylab("Proportion of constraint") +
  #scale_color_manual(values=c("grey","orange")) +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line=element_line(colour = "black"))
enhancer_constraint_large_intergenic_regions_plot

td_enhancer_set <- enhancers_large_intergenic_regions$TD_associated == "TRUE"
non_td_enhancer_set <- enhancers_large_intergenic_regions$TD_associated == "FALSE"
enhancer_td_mean <- mean(enhancers_large_intergenic_regions$prop_constraint[td_enhancer_set])
enhancer_non_td_mean <- mean(enhancers_large_intergenic_regions$prop_constraint[non_td_enhancer_set])

enhancer_constraint_means <- data.frame(row.names = c(1,2))
enhancer_constraint_means$type = c("trans-dev genes", "other")
enhancer_constraint_means$mean = c(enhancer_td_mean, enhancer_non_td_mean)


# Calculate 95% confidence intervals using bootstrapping

enhancer_boot_out = list() # Create lists for storing bootstrap outputs
enhancer_boot_ci_out = list()

enhancer_cats=c("TRUE","FALSE")

# Loop over each category and calculate bootstrap means and CIs of constraint

for (cat in 1:length(enhancer_cats)) {
  
  temp=enhancers_large_intergenic_regions[enhancers_large_intergenic_regions$TD_associated==enhancer_cats[cat],]
  enhancer_boot_out[[cat]] <- boot(temp$prop_constraint, samplemean, R=1000)
  enhancer_boot_ci_out[[cat]] <- boot.ci(enhancer_boot_out[[cat]], type = "perc")
  # Extract confidence intervals and add to kin_means df
  enhancer_constraint_means$lower.ci[cat] <- (enhancer_boot_ci_out[[cat]][["percent"]][(4)])
  enhancer_constraint_means$upper.ci[cat] <- (enhancer_boot_ci_out[[cat]][["percent"]][(5)])

}

# Set the order of categories for plotting
enhancer_constraint_means$type = factor(enhancer_constraint_means$type, levels = c("trans-dev genes", "other"))
# And plot
enhancer_mean_constraint <- ggplot(data=enhancer_constraint_means, aes(x = type, y = mean)) + 
  geom_col(width = 0.5, fill="lightgrey") +
  geom_errorbar(ymin=enhancer_constraint_means$lower.ci, ymax=enhancer_constraint_means$upper.ci, width=0.2) +
  scale_y_continuous(expand=c(0,0), limits=c(0,0.125)) +
  xlab("Developmental transcription factor neighbour?") +
  ylab("Proportion of constraint") +
  scale_x_discrete(labels=c("Yes", "No")) +
  scale_color_manual(values=c("grey","orange")) +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line=element_line(colour = "black"))
enhancer_mean_constraint

# Test for significant difference
t.test(enhancers_large_intergenic_regions$prop_constraint[td_enhancer_set],enhancers_large_intergenic_regions$prop_constraint[non_td_enhancer_set])
wilcox.test(enhancers_large_intergenic_regions$prop_constraint[td_enhancer_set],enhancers_large_intergenic_regions$prop_constraint[non_td_enhancer_set])
length(enhancers_large_intergenic_regions$prop_constraint[td_enhancer_set])
length(enhancers_large_intergenic_regions$prop_constraint[non_td_enhancer_set])

# Plot intergenic and enhancer constraint plots together
noncoding_constraint_plots <- grid.arrange(intergen_top5_mean_constraint, enhancer_mean_constraint,
                                                 nrow = 2, ncol=1)
# Add labels to the arranged plots
noncoding_constraint_plots_labelled <- as_ggplot(noncoding_constraint_plots) +                                # transform to a ggplot
  draw_plot_label(label = c("A", "B"), size = 15,
                  x = c(0,0), y = c(1,0.5)) # Add labels
noncoding_constraint_plots_labelled 


## Look at density of td genes per chromosome and compare to overall chromosome constraint
# Import trans-dev gene coords
td_genes <- read_delim("trans_dev_gene_coords.bed","\t", escape_double = FALSE, trim_ws = TRUE)
td_genes$length <- td_genes$stop - td_genes$start


chroms = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11",
           "chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22")

td_gene_lengths <- data.frame(row.names = c(1:22))
td_gene_lengths$chromosome = chroms

for (chrom in 1:length(chroms)) {
  temp=td_genes[td_genes$chr==chroms[chrom],] # This will create a subset for each chromosome
  td_gene_lengths$total_td_length[chrom] <- sum(temp$length)
  td_gene_lengths$total_genes[chrom] <- length(temp$chr)
}

#Add in chr lengths  
chr_lengths <- read_delim("hg38_chr_lengths.txt","\t", escape_double = FALSE, trim_ws = TRUE)
td_gene_lengths$chr_length <- chr_lengths$length
td_gene_lengths$prop_td <- td_gene_lengths$total_td_length/td_gene_lengths$chr_length
td_gene_lengths$bp_per_td <- td_gene_lengths$chr_length/td_gene_lengths$total_genes

# Compare to constraint on chromosome
chr_constraint <- read_delim("chromosome_constraint.txt","\t", escape_double = FALSE, trim_ws = TRUE)
td_gene_lengths$chr_constraint <- chr_constraint$prop_constraint

td_gene_lengths$chromosome <- factor(td_gene_lengths$chromosome, levels =  c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11",
                               "chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22"))

# Plots

# Number of td genes on chromosome

cor.test(td_gene_lengths$total_genes, td_gene_lengths$chr_constraint)

total_genes_plot <- ggplot(data=td_gene_lengths, aes(x = total_genes, y = chr_constraint, colour=chromosome)) + 
  geom_point() +
  annotate("text", label = "Pearson's r = 0.47, p = 0.03", x = 200, y = 0.24, size = 4, colour = "grey3") +
  scale_y_continuous(expand=c(0,0), limits=c(0,0.25)) +
  xlab("Number of trans-dev genes on chromosome") +
  ylab("Chromosome constraint") +
  #scale_color_manual(values=c("grey","orange")) +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line=element_line(colour = "black"), legend.position="BLANK")
total_genes_plot



# Proportion of chromosome that is td genes

cor.test(td_gene_lengths$prop_td, td_gene_lengths$chr_constraint)

prop_td_plot <- ggplot(data=td_gene_lengths, aes(x = prop_td, y = chr_constraint, colour=chromosome)) + 
  geom_point() +
  annotate("text", label = "Pearson's r = 0.07, p = 0.8", x = 0.08, y = 0.22, size = 4, colour = "grey3") +
  scale_y_continuous(expand=c(0,0), limits=c(0,0.25)) +
  xlab("Proportion of chromosome that is trans-dev genes") +
  ylab("Chromosome constraint") +
  guides(colour=guide_legend(ncol=3,title="Chromosome")) +
  #scale_color_manual(values=c("grey","orange")) +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line=element_line(colour = "black"),
        legend.box.background = element_rect(colour = "black"))
prop_td_plot

# base pairs per td

cor.test(td_gene_lengths$bp_per_td, td_gene_lengths$chr_constraint)

bp_per_td_plot <- ggplot(data=td_gene_lengths, aes(x = bp_per_td/1000000, y = chr_constraint, colour=chromosome)) + 
  geom_point() +
  scale_y_continuous(expand=c(0,0), limits=c(0,0.25)) +
  annotate("text", label = "Pearson's r = -0.46, p = 0.03", x = 2, y = 0.22, size = 4, colour = "grey3") +
  xlab("Base pairs per trans-dev gene (Mbp)") +
  ylab("Chromosome constraint") +
  #scale_color_manual(values=c("grey","orange")) +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line=element_line(colour = "black"), legend.position="BLANK")
bp_per_td_plot

td_constraint_figs <- grid.arrange(prop_td_plot,bp_per_td_plot,total_genes_plot,
                     nrow = 2, ncol=2,layout_matrix = rbind(c(1,1), c(2,3)))
# # Add labels to the arranged plots
# fig4_labelled <- as_ggplot(fig4 ) +                                # transform to a ggplot
#   draw_plot_label(label = c("A", "B", "C","D","E"), size = 15,
#                   x = c(0,0.33,0.66,0, 0.6), y = c(1, 1, 1,0.5,0.5)) # Add labels
# fig4_labelled              

# But does this just reflect number of genes, trans-dev or not, explains constraint?

genes_per_chr <- read_delim("number_of_genes_per_chr.txt","\t", escape_double = FALSE, trim_ws = TRUE)
genes_per_chr$constraint <- chr_constraint$prop_constraint
genes_per_chr$td_genes <- td_gene_lengths$total_genes

cor.test(genes_per_chr$n.genes, genes_per_chr$constraint)

genes_per_chr_plot <- ggplot(data=genes_per_chr, aes(x = n.genes, y = constraint, colour=chr)) + 
  geom_point() +
  #annotate("text", label = "Pearson's r = 0.47, p = 0.03", x = 200, y = 0.24, size = 4, colour = "grey3") +
  #scale_y_continuous(expand=c(0,0), limits=c(0,0.25)) +
  xlab("Number of genes on chromosome") +
  ylab("Chromosome constraint") +
  #scale_color_manual(values=c("grey","orange")) +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line=element_line(colour = "black"))
genes_per_chr_plot


cor.test(genes_per_chr$n.genes, genes_per_chr$td_genes)

genes_per_chr_v_td_genes_plot <- ggplot(data=genes_per_chr, aes(x = n.genes, y = td_genes, colour=chr)) + 
  geom_point() +
  #annotate("text", label = "Pearson's r = 0.47, p = 0.03", x = 200, y = 0.24, size = 4, colour = "grey3") +
  #scale_y_continuous(expand=c(0,0), limits=c(0,0.25)) +
  xlab("Number of genes on chromosome") +
  ylab("Number of trans-dev genes") +
  #scale_color_manual(values=c("grey","orange")) +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line=element_line(colour = "black"))
genes_per_chr_v_td_genes_plot

