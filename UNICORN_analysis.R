library(tidyverse)
library(ggpubr)
library(gridExtra)
library(cowplot)
library(ggbiplot)
library(unix)
library(GGally)
library(dplyr)
library(ggpointdensity)
library(viridis)

options(scipen=10000)

# This script is for analysis of UNannotated Intergenic COnstraint RegioNs (UNICORNs)


UNICORNs <- read_delim("UNICORNs_FINAL_set.bed", 
                         "\t", escape_double = FALSE, col_names = TRUE, 
                         trim_ws = TRUE)
UNICORNs$length <- UNICORNs$end-UNICORNs$start
UNICORN_mean_length <- mean(UNICORNs$length)
UNICORN_sd <- sd(UNICORNs$length)

quantile(UNICORNs$length, probs=seq(0,1,0.95))

ggplot(data=UNICORNs, aes(x=length)) +
  geom_histogram(binwidth = 10) +
  xlab("Size of UNICORNs") +
  #geom_vline(xintercept = 482.881, lty=3,colour="red") +
  theme_bw() + 
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16,face="bold")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))

ggplot(data=UNICORNs, aes(x=region,y=length)) +
  #geom_violin() +
  geom_boxplot() +
  theme_bw() + 
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16,face="bold")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))

# Zoom in to tail of plot
ggplot(data=UNICORNs, aes(x=length)) +
  geom_histogram(binwidth = 10) +
  xlab("Size of UNICORN") +
  #geom_vline(xintercept = 482.881, lty=3,colour="red") +
  scale_x_continuous(limits=c(200,1500)) +
  theme_bw() + 
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16,face="bold")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))

# Take the top 99.9th percentile and import into GREAT
quantile(UNICORN_master_df$length, probs=seq(0,1,0.999))
UNICORN_99.9th_pc_length <- subset(UNICORN_master_df, length>=514)

UNICORN_99.9th_pc_coords <- data.frame(UNICORN_99.9th_pc_length$intergen_chrom,UNICORN_99.9th_pc_length$start,UNICORN_99.9th_pc_length$end)
write.table(UNICORN_99.9th_pc_coords, file="UNICORN_99.9th_pc_coords.bed", sep="\t", quote=FALSE, row.names = FALSE)




# Test for enrichment of trans-dev genes in the gene set from GREAT associated with the long intergenic regions of
# high constraint:
# Create 2x2 matrix
intergenic_contingency_table <- rbind(c(19293,2863),c(378,100))
colnames(intergenic_contingency_table) <- c("all_genes","td_genes")
rownames(intergenic_contingency_table) <- c("whole_genome","intergenic_constraint")
fisher.test(intergenic_contingency_table)
# Fisher's Exact Test for Count Data
# 
# data:  intergenic_contingency_table
# p-value = 0.000001401
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  1.410504 2.236239
# sample estimates:
# odds ratio 
#   1.782676 


# Identify nearest genes to each element using bedtools closest:
# bedtools closest -a intergenic_constraint_regions_gt_10bp_position_defined_sorted.bed -b gencode.v37.exons.sorted.bed -d > intergenic_constraint_regions_gt_10bp_nearest_genes.txt


# Bring in features: GWAS catolog SNPs, TOPMed SNP counts and AFs, recombination rate, distance to regulatory features and genes
# Import back in and look at relationship between size of elements and proximity to genes

UNICORN_exon_proximity <- read_delim("./distance_to_annotated_features/UNICORNs_closest_exons.bed", 
                                      "\t", escape_double = FALSE, col_names = TRUE, 
                                      trim_ws = TRUE)

dhs <- read_delim("./distance_to_annotated_features/UNICORNs_closest_dhs.bed", 
                  "\t", escape_double = FALSE, trim_ws = TRUE)
ccre <- read_delim("./distance_to_annotated_features/UNICORNs_closest_cCREs.bed", 
                   "\t", escape_double = FALSE, trim_ws = TRUE)
anchor <- read_delim("./distance_to_annotated_features/UNICORNs_closest_chromatin_anchor.bed", 
                     "\t", escape_double = FALSE, trim_ws = TRUE)

GWAS_SNP_counts <- read_delim("./GWAS_signals/UNICORN_GWAS_SNP_counts.bed", 
                              "\t", escape_double = FALSE, trim_ws = TRUE)

TOPMed_SNPs <- read_delim("./SNPs/UNICORN_SNPs.bed", 
                          "\t", escape_double = FALSE, trim_ws = TRUE)

#recombination <- read_delim("./analysis/Recombination/intergenic_constraint_regions_gt_10bp_recombination_rates.txt", 
#                           "\t", escape_double = FALSE, trim_ws = TRUE)

UNICORN_master_df <- merge(x=UNICORN_exon_proximity[,c("U_chrom","U_start","U_end","length","type","ID","name","distance")],y=GWAS_SNP_counts[,c("U_chrom","U_start","GWAS_SNPs_count")],by=c("U_chrom","U_start"))
# Remove duplicate rows (not sure why these exist...)
UNICORN_master_df <- UNICORN_master_df %>% distinct()

UNICORN_master_df <- merge(x=UNICORN_master_df,y=dhs[,c("U_chrom","U_start","dhs_distance")],by=c("U_chrom","U_start"))
UNICORN_master_df <- merge(x=UNICORN_master_df,y=ccre[,c("U_chrom","U_start","cCRE_distance")],by=c("U_chrom","U_start"))
UNICORN_master_df <- merge(x=UNICORN_master_df,y=anchor[,c("U_chrom","U_start","anchor_distance")],by=c("U_chrom","U_start"))
UNICORN_master_df <- merge(x=UNICORN_master_df,y=TOPMed_SNPs[,c("U_chrom","U_start","av_MAF","SNPs_per_bp")],by=c("U_chrom","U_start"))
#UNICORN_master_df <- merge(x=UNICORN_master_df,y=recombination[,c("intergen_chrom","start","recombination_rate")],by=c("intergen_chrom","start"))

# Rename columns to match those of nonconstraint df
names(UNICORN_master_df)[names(UNICORN_master_df) == "U_chrom"] <- "intergen_chrom"
names(UNICORN_master_df)[names(UNICORN_master_df) == "U_start"] <- "start"
names(UNICORN_master_df)[names(UNICORN_master_df) == "U_end"] <- "end"
names(UNICORN_master_df)[names(UNICORN_master_df) == "distance"] <- "gene_distance"
names(UNICORN_master_df)[names(UNICORN_master_df) == "ID"] <- "gene_ID"
names(UNICORN_master_df)[names(UNICORN_master_df) == "name"] <- "gene_name"
names(UNICORN_master_df)[names(UNICORN_master_df) == "cCRE_distance"] <- "ccre_distance"

# Look at how the different features relate:

ggpairs(UNICORN_master_df[,c(4,7:14)])


# Output UNICORN dataset
write.table(UNICORN_master_df, file="UNICORN_data.txt", sep="\t",quote=FALSE, row.names=FALSE)


# Bring in datasets for non-constraint regions
setwd("~/Google Drive/UU_postdoc/KLT_Lab/Zoonomia/unannotated_constraint/intergenic_regions/analysis/non_constraint_regions/distance_to_annotated_features")
non_constraint_gene_distance <- read_delim("./intergenic_nonconstraint_clusters_nearest_genes.bed", "\t", escape_double = FALSE, trim_ws = TRUE)
non_constraint_gene_distance <- subset(non_constraint_gene_distance, gene_distance>1000)
non_constraint_ccre_distance <- read_delim("./intergenic_nonconstraint_clusters_nearest_ccre.bed", "\t", escape_double = FALSE, trim_ws = TRUE)
non_constraint_dhs_distance <- read_delim("./intergenic_nonconstraint_clusters_nearest_dhs.bed", "\t", escape_double = FALSE, trim_ws = TRUE)
non_constraint_anchor_distance <- read_delim("./intergenic_nonconstraint_clusters_nearest_chromatin_anchor.bed", "\t", escape_double = FALSE, trim_ws = TRUE)

setwd("~/Google Drive/UU_postdoc/KLT_Lab/Zoonomia/unannotated_constraint/intergenic_regions/analysis/non_constraint_regions/GWAS_signals")
non_constraint_GWAS <- read_delim("./intergenic_nonconstraint_clusters_GWAS_SNP_counts.bed", "\t", escape_double = FALSE, trim_ws = TRUE)
setwd("~/Google Drive/UU_postdoc/KLT_Lab/Zoonomia/unannotated_constraint/intergenic_regions/analysis/non_constraint_regions/SNPs")
non_constraint_SNPs <- read_delim("./intergenic_nonconstraint_clusters_SNPs.bed", "\t", escape_double = FALSE, trim_ws = TRUE)

setwd("~/Google Drive/UU_postdoc/KLT_Lab/Zoonomia/unannotated_constraint/intergenic_regions/analysis/non_constraint_regions/recombination")
non_constraint_recombination <-read_delim("./intergenic_nonconstraint_clusters_recombination_rates.bed", "\t", escape_double = FALSE, trim_ws = TRUE)

nonUNICORN_master_df <- merge(x=non_constraint_gene_distance[,c("intergen_chrom","start","end","length","gene_name","gene_distance")],y=non_constraint_ccre_distance[,c("intergen_chrom","start","ccre_distance")],by=c("intergen_chrom","start"))
nonUNICORN_master_df <- merge(x=nonUNICORN_master_df,y=non_constraint_dhs_distance[,c("intergen_chrom","start","dhs_distance")],by=c("intergen_chrom","start"))
nonUNICORN_master_df <- merge(x=nonUNICORN_master_df,y=non_constraint_anchor_distance[,c("intergen_chrom","start","anchor_distance")],by=c("intergen_chrom","start"))
nonUNICORN_master_df <- merge(x=nonUNICORN_master_df,y=non_constraint_GWAS[,c("intergen_chrom","start","GWAS_SNPs_count")],by=c("intergen_chrom","start"))
nonUNICORN_master_df <- merge(x=nonUNICORN_master_df,y=non_constraint_SNPs[,c("intergen_chrom","start","av_MAF","SNPs_per_bp")],by=c("intergen_chrom","start"))
nonUNICORN_master_df <- merge(x=nonUNICORN_master_df,y=non_constraint_recombination[,c("intergen_chrom","start","recombination_rate")],by=c("intergen_chrom","start"))


# Remove any row containing 'NA'
nonUNICORN_master_df <- na.omit(nonUNICORN_master_df)

# Set max length of clusters to be the same as constraint clusters (1325)
nonUNICORN_master_df_subset <- subset(nonUNICORN_master_df, length<1325)

# Randomly sample the df to reduce its size to the same as the constraint set
nonUNICORN_master_df_subset <- sample_n(nonUNICORN_master_df_subset, 424179)

write.table(nonUNICORN_master_df_subset, file="unannotated_intergenic_regions_non_constraint.bed", sep="\t",row.names=FALSE, quote=FALSE)


# combine dfs after adding definition column to both
UNICORN_master_df$constraint_type = "constraint"
nonUNICORN_master_df_subset$constraint_type = "non_constraint"

# Make sure both dfs have same columns and are labelled the same, then combine for plotting:

UNICORN_master_df_13_vars <- subset(UNICORN_master_df, select = -c(5,6))
nonUNICORN_master_df_13_vars <- subset(nonUNICORN_master_df_subset, select = -c(13))

cluster_features_master <- rbind(UNICORN_master_df_13_vars,nonUNICORN_master_df_13_vars)


# Compare size of clusters with a histogram

cluster_size_histogram <- ggplot(data=cluster_features_master, aes(x=length, fill=constraint_type)) +
  geom_histogram(binwidth = 10,alpha=0.5, position="identity",aes(y = (..count..)/sum(..count..))) +
  xlab("Size of unannotated intergenic cluster") +
  ylab("Frequency") +
  labs(fill="Constraint category") +
  #geom_vline(xintercept = 482.881, lty=3,colour="red") +
  #scale_x_continuous(limits=c(0,1500)) +
  scale_fill_manual(values=c("firebrick3","grey1"), labels=c("Constraint","Non-constraint")) +
  theme_bw() + 
  theme(axis.text=element_text(size=10),axis.title=element_text(size=12,face="bold")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),legend.position=c(0.6,0.6),axis.title.x = element_blank())
cluster_size_histogram
# Doesn't look too bad, quite a lot more smaller clusters in constraint set though


# Plot AF v SNPs_per_bp
cluster_features_master %>%
  subset(length>200) %>%
  subset(constraint_type=="constraint") %>%
  subset(constraint_type=="non_constraint") %>%
  arrange(desc(constraint_type)) %>%
ggplot(aes(x=log10(av_MAF),y=SNPs_per_bp, colour=constraint_type, alpha=constraint_type)) +
  geom_point(size=0.6) +
  #geom_pointdensity() +
  #scale_color_viridis() +
  xlab("Average minor allele frequency") +
  ylab("SNPs per bp") +
  #labs(colour="") +
  #geom_vline(xintercept = 482.881, lty=3,colour="red") +
  #scale_x_continuous(limits=c(0,1500)) +
  scale_alpha_manual(guide='none', values = c(non_constraint = 1, constraint = 0.2)) +
  scale_colour_manual(values=c("firebrick3","darkgrey"), labels=c("UNICORNs","Background")) +
  theme_bw() + 
  theme(axis.text=element_text(size=10),axis.title=element_text(size=12,face="bold")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),legend.position=c(0.8,0.2))

# Show percentiles on plots
constraint_DAF_mean <- cluster_features_master %>%
  subset(length>200) %>%
  subset(constraint_type=="constraint") %>%
  summarise(mean(av_MAF))
constraint_DAF_sd <- cluster_features_master %>%
  subset(length>200) %>%
  subset(constraint_type=="constraint") %>%
  summarise(sd(av_MAF))

cluster_features_master %>%
  subset(length>200) %>%
  subset(constraint_type=="constraint") %>%
  summarise(x = quantile(av_MAF, c(0.25, 0.75)), q = c(0.25, 0.75))

# point density
constraint_dotplot <- cluster_features_master %>%
  subset(length>200) %>%
  subset(constraint_type=="constraint") %>%
  #subset(constraint_type=="non_constraint") %>%
  arrange(desc(constraint_type)) %>%
  ggplot(aes(x=log10(av_MAF),y=SNPs_per_bp)) +
  #geom_point(size=0.6) +
  geom_pointdensity(size=0.2) +
  scale_color_viridis() +
  #geom_vline(xintercept=log10(0.00001596208), lty=2, colour="darkgrey") + #bottom 25%
  #geom_vline(xintercept=log10(0.00050802427), lty=2, colour="darkgrey") + #top 75%
  xlab("Average minor allele frequency") +
  ylab("") +
  #labs(colour="") +
  #geom_vline(xintercept = 482.881, lty=3,colour="red") +
  scale_x_continuous(limits=c(-8,0), labels=c("10^-8",10^-6,10^-4,10^-2,10^0)) +
  scale_y_continuous(limits=c(0,0.8)) +
  theme_bw() + 
  theme(axis.text=element_text(size=10),axis.title=element_text(size=12,face="bold")) +
  theme(panel.border = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),legend.position=c(0.1,0.7), axis.line.x = element_blank(),
        axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), legend.title = element_blank(),)

cluster_features_master %>%
  subset(length>200) %>%
  subset(constraint_type=="non_constraint") %>%
  summarise(x = quantile(av_MAF, c(0.25, 0.75)), q = c(0.25, 0.75))

non_constraint_dotplot <- cluster_features_master %>%
  subset(length>200) %>%
  #subset(constraint_type=="constraint") %>%
  subset(constraint_type=="non_constraint") %>%
  arrange(desc(constraint_type)) %>%
  ggplot(aes(x=log10(av_MAF),y=SNPs_per_bp)) +
  #geom_point(size=0.6) +
  geom_pointdensity(size=0.2) +
  scale_color_viridis() +
  #geom_vline(xintercept=log10(0.00006942993), lty=2, colour="darkgrey") + #bottom 25%
  #geom_vline(xintercept=log10(0.00146638518), lty=2, colour="darkgrey") + #top 75%
  xlab("Average derived allele frequency") +
  ylab("") +
  #labs(colour="") +
  #geom_vline(xintercept = 482.881, lty=3,colour="red") +
  scale_x_continuous(limits=c(-8,0), labels=c(expression(10^-8),expression(10^-6),expression(10^-4),expression(10^-2),1)) +
  scale_y_continuous(limits=c(0,0.8)) +
  theme_bw() + 
  theme(axis.text=element_text(size=10),axis.title=element_text(size=12,face="bold")) +
  theme(panel.border = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),legend.position="none")

dotplots <- ggarrange(constraint_dotplot,non_constraint_dotplot,nrow=2,ncol=1, heights=c(1,1.2))

# Add common y-axis
annotate_figure(dotplots, left = text_grob("SNPs per bp", color = "black", rot = 90, 
                                           face = "bold", size = 12))

# Bin the data for plotting
summary(cluster_features_master)
# set up cut-off values 
breaks <- c(10,20,50,100,500,1325)
# specify interval/bin labels
tags <- c("11-20","21-50", "51-100", "101-500", "501-1325")
# bucketing values into bins
group_tags <- cut(cluster_features_master$length, 
                  breaks=breaks, 
                  include.lowest=TRUE, 
                  right=TRUE, 
                  labels=tags)
# inspect bins
summary(group_tags)

cluster_features_master_binned <- as_tibble(cluster_features_master) %>% 
  mutate(tag = case_when(
    length <= 20 ~ tags[1],
    length >= 21 & length <= 50 ~ tags[2],
    length >= 51 & length <= 100 ~ tags[3],
    length >= 101 & length <= 500 ~ tags[4],
    length >= 501 & length < 2000 ~ tags[5]
  ))

# Make tag a factor
cluster_features_master_binned$tag <- factor(cluster_features_master_binned$tag,
                                             levels = tags,
                                             ordered = FALSE)

# Set constraint and nonconstraint entries as factors for running stats
constraint<- cluster_features_master_binned$constraint_type=="constraint"
nonconstraint<- cluster_features_master_binned$constraint_type=="non_constraint"

cluster_features_master_binned$constraint_binomial = 0
cluster_features_master_binned$constraint_binomial[constraint] = 1


# Boxplots
# Distance to gene
gene_distance_plot <- ggplot(data = cluster_features_master_binned, aes(x=tag,y=gene_distance/1000,colour=constraint_type, fill=constraint_type)) + 
  #geom_jitter(alpha=0.2) +
  geom_boxplot(alpha=0.3,outlier.shape=NA) + 
  labs(x='Intergenic cluster length', y='Distance to gene (kbp)') +
  #guides(color=FALSE) +
  scale_colour_manual(values=c("tomato","grey2")) +
  scale_fill_manual(values=c("firebrick3","grey1")) +
  theme_bw() + 
  theme(axis.text=element_text(size=10),axis.title=element_text(size=12,face="bold")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),legend.position = "none",axis.title.x = element_blank(),axis.text.x = element_blank()) +
  coord_cartesian(ylim=c(0,150))

gene_distance_plot

t.test(cluster_features_master_binned$gene_distance[constraint], cluster_features_master_binned$gene_distance[nonconstraint])
lapply(split(cluster_features_master_binned, cluster_features_master_binned$tag), function(i){plot(lm(constraint_type~gene_distance, data=i))})
lapply(split(cluster_features_master_binned, cluster_features_master_binned$tag), function(i){anova(lm(constraint_binomial~gene_distance, data=i))})


# Distance to cCRE

ccre_distance_plot <- ggplot(data = cluster_features_master_binned, aes(x=tag,y=ccre_distance/1000,colour=constraint_type, fill=constraint_type)) + 
  #geom_jitter(alpha=0.2) +
  geom_boxplot(alpha=0.3,outlier.shape=NA) + 
  labs(x='Intergenic cluster length', y='Distance to  cCRE (kbp)') +
  #guides(color=FALSE) +
  scale_colour_manual(values=c("tomato","grey2")) +
  scale_fill_manual(values=c("firebrick3","grey1")) +
  theme_bw() + 
  theme(axis.text=element_text(size=10),axis.title=element_text(size=12,face="bold")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_blank()) +
  coord_cartesian(ylim=c(0,16))

ccre_distance_plot

t.test(cluster_features_master_binned$ccre_distance[constraint], cluster_features_master_binned$ccre_distance[nonconstraint])
lapply(split(cluster_features_master_binned, cluster_features_master_binned$tag), function(i){plot(lm(constraint_type~ccre_distance, data=i))})
lapply(split(cluster_features_master_binned, cluster_features_master_binned$tag), function(i){anova(lm(constraint_binomial~ccre_distance, data=i))})


# Distance to DHS

dhs_distance_plot <- ggplot(data = cluster_features_master_binned, aes(x=tag,y=dhs_distance/1000,colour=constraint_type, fill=constraint_type)) + 
  #geom_jitter(alpha=0.2) +
  geom_boxplot(alpha=0.3,outlier.shape=NA) + 
  labs(x='Intergenic cluster length (bp)', y='Distance to DHS (kb)') +
  #guides(color=FALSE) +
  scale_colour_manual(values=c("tomato","grey2")) +
  scale_fill_manual(values=c("firebrick3","grey1")) +
  theme_bw() + 
  theme(axis.text=element_text(size=10),axis.title=element_text(size=12,face="bold")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),legend.position = "none") +
  coord_cartesian(ylim=c(0,6))

dhs_distance_plot

t.test(cluster_features_master_binned$dhs_distance[constraint], cluster_features_master_binned$dhs_distance[nonconstraint])
lapply(split(cluster_features_master_binned, cluster_features_master_binned$tag), function(i){lm(constraint_binomial~dhs_distance, data=i)})
lapply(split(cluster_features_master_binned, cluster_features_master_binned$tag), function(i){anova(lm(constraint_binomial~dhs_distance, data=i))})

bin_51_100_dhs <- subset(cluster_features_master_binned, tag=="51-100")
dhs_constraint <- bin_51_100_dhs$constraint_type == "constraint"
dhs_nonconstraint <- bin_51_100_dhs$constraint_type == "non_constraint"
mean(bin_51_100_dhs $dhs_distance[dhs_constraint])
mean(bin_51_100_dhs $dhs_distance[dhs_nonconstraint])

mean(bin_11_20$dhs_distance[bin_11_20_constraint])
mean(bin_11_20$dhs_distance[bin_11_20_nonconstraint])
mean(bin_21_50$dhs_distance[bin_21_50_constraint])
mean(bin_21_50$dhs_distance[bin_21_50_nonconstraint])
mean(bin_51_100$dhs_distance[bin_51_100_constraint])
mean(bin_51_100$dhs_distance[bin_51_100_nonconstraint])
mean(bin_101_500$dhs_distance[bin_101_500_constraint])
mean(bin_101_500$dhs_distance[bin_101_500_nonconstraint])
mean(bin_501_1325$dhs_distance[bin_501_1325_constraint])
mean(bin_501_1325$dhs_distance[bin_501_1325_nonconstraint])

#Distance to chromatin loop anchor

loop_anchor_distance_plot <- ggplot(data = cluster_features_master_binned, aes(x=tag,y=anchor_distance/1000,colour=constraint_type, fill=constraint_type)) + 
  #geom_jitter(alpha=0.2) +
  geom_boxplot(alpha=0.3,outlier.shape=NA) + 
  labs(x='Intergenic cluster length', y='Distance to loop anchor (Kbp)') +
  #guides(color=FALSE) +
  scale_colour_manual(values=c("tomato","grey2")) +
  scale_fill_manual(values=c("firebrick3","grey1")) +
  theme_bw() + 
  theme(axis.text=element_text(size=10),axis.title=element_text(size=12,face="bold")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),legend.position = "none", axis.title.x = element_blank(),axis.text.x = element_blank()) +
  coord_cartesian(ylim=c(0,120))

loop_anchor_distance_plot

t.test(cluster_features_master_binned$anchor_distance[constraint], cluster_features_master_binned$anchor_distance[nonconstraint])
lapply(split(cluster_features_master_binned, cluster_features_master_binned$tag), function(i){plot(lm(constraint_type~anchor_distance, data=i))})
lapply(split(cluster_features_master_binned, cluster_features_master_binned$tag), function(i){anova(lm(constraint_binomial~anchor_distance, data=i))})


# GWAS SNPs count
GWAS_SNPs_plot <- ggplot(data = cluster_features_master_binned, aes(x=tag,y=GWAS_SNPs_count,colour=constraint_type, fill=constraint_type)) + 
  #geom_jitter(alpha=0.2) +
  geom_boxplot(alpha=0.3) + 
  labs(x='Intergenic cluster length', y='GWAS hits') +
  #guides(color=FALSE) +
  scale_colour_manual(values=c("tomato","grey2")) +
  scale_fill_manual(values=c("firebrick3","grey1")) +
  theme_bw() + 
  theme(axis.text=element_text(size=10),axis.title=element_text(size=12,face="bold")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),legend.position = "none",axis.title.x = element_blank(),axis.text.x = element_blank()) #+
#coord_cartesian(ylim=c(0,10))

GWAS_SNPs_plot

# t.test(cluster_features_master_binned$GWAS_SNPs_count[constraint], cluster_features_master_binned$GWAS_SNPs_count[nonconstraint])
# lapply(split(cluster_features_master_binned, cluster_features_master_binned$tag), function(i){plot(lm(constraint_type~GWAS_SNPs_count, data=i))})
# lapply(split(cluster_features_master_binned, cluster_features_master_binned$tag), function(i){anova(lm(constraint_binomial~GWAS_SNPs_count, data=i))})

# Wilcoxon rank sum test
#Overall
wilcox.test(GWAS_SNPs_count~constraint_type, data=cluster_features_master,alternative = "less")

# And per bin
lapply(split(cluster_features_master_binned, cluster_features_master_binned$tag), function(i){wilcox.test(GWAS_SNPs_count~constraint_type, data=i,alternative = "less")})
# Bonferonni correction
#bin 11-20
p.adjust(0.000000000005533, method = "bonferroni",n = 5) # pcorr = 0.000000000027665
#bin 21-50
p.adjust(0.0000000001282, method = "bonferroni",n = 5) # p-corr = 0.000000000641
#bin 51-100
p.adjust(0.000003308, method = "bonferroni",n = 5) # p-corr = 0.00001654
#bin 101-500
p.adjust(0.00000000000000022, method = "bonferroni",n = 5) # p-corr = 0.0000000000000011
#bin 501-1325
p.adjust(0.09536, method = "bonferroni",n = 5) # p-corr = 0.4768

# Average MAF
av_MAF_plot <- ggplot(data = cluster_features_master_binned, aes(x=tag,y=av_MAF,colour=constraint_type, fill=constraint_type)) + 
  #geom_jitter(alpha=0.2) +
  geom_boxplot(alpha=0.3,outlier.shape=NA) + 
  labs(x='Intergenic cluster length (bp)', y='MAF') +
  #guides(color=FALSE) +
  scale_colour_manual(values=c("tomato","grey2")) +
  scale_fill_manual(values=c("firebrick3","grey1")) +
  theme_bw() + 
  theme(axis.text=element_text(size=10),axis.text.x=element_text(angle = 45,hjust = 1),axis.title=element_text(size=12,face="bold")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),legend.position = "none") +
  coord_cartesian(ylim=c(0,0.005))

av_MAF_plot

# Test for significant differences:
# Wilcoxon rank sum test
#Overall
wilcox.test(av_MAF~constraint_type, data=cluster_features_master,alternative = "less")

# And per bin
lapply(split(cluster_features_master_binned, cluster_features_master_binned$tag), function(i){wilcox.test(av_MAF~constraint_type, data=i,alternative = "less")})

# Bonferonni correction
#All bins have p = 0.00000000000000022, so same corrected p-value for each
p.adjust(0.00000000000000022, method = "bonferroni",n = 5) # pcorr = 0.0000000000000011


t.test(cluster_features_master_binned$av_MAF[constraint], cluster_features_master_binned$av_MAF[nonconstraint])
lapply(split(cluster_features_master_binned, cluster_features_master_binned$tag), function(i){plot(lm(constraint_type~av_MAF, data=i))})
lapply(split(cluster_features_master_binned, cluster_features_master_binned$tag), function(i){anova(lm(constraint_binomial~av_MAF, data=i))})

bin_11_20 <- subset(cluster_features_master_binned, tag=="11-20")
bin_11_20_constraint <- bin_11_20$constraint_type == "constraint"
bin_11_20_nonconstraint <- bin_11_20$constraint_type == "non_constraint"
mean(bin_11_20$av_MAF[bin_11_20_constraint])
mean(bin_11_20$av_MAF[bin_11_20_nonconstraint])

bin_21_50 <- subset(cluster_features_master_binned, tag=="21-50")
bin_21_50_constraint <- bin_21_50$constraint_type == "constraint"
bin_21_50_nonconstraint <- bin_21_50$constraint_type == "non_constraint"
mean(bin_21_50$av_MAF[bin_21_50_constraint])
mean(bin_21_50$av_MAF[bin_21_50_nonconstraint])

bin_51_100 <- subset(cluster_features_master_binned, tag=="51-100")
bin_51_100_constraint <- bin_51_100$constraint_type == "constraint"
bin_51_100_nonconstraint <- bin_51_100$constraint_type == "non_constraint"
mean(bin_51_100$av_MAF[bin_51_100_constraint])
mean(bin_51_100$av_MAF[bin_51_100_nonconstraint])

bin_101_500 <- subset(cluster_features_master_binned, tag=="101-500")
bin_101_500_constraint <- bin_101_500$constraint_type == "constraint"
bin_101_500_nonconstraint <- bin_101_500$constraint_type == "non_constraint"
mean(bin_101_500$av_MAF[bin_101_500_constraint])
mean(bin_101_500$av_MAF[bin_101_500_nonconstraint])


bin_501_1325 <- subset(cluster_features_master_binned, tag=="501-1325")
bin_501_1325_constraint <- bin_501_1325$constraint_type == "constraint"
bin_501_1325_nonconstraint <- bin_501_1325$constraint_type == "non_constraint"
mean(bin_501_1325$av_MAF[bin_501_1325_constraint])
mean(bin_501_1325$av_MAF[bin_501_1325_nonconstraint])


# SNPs per bp
SNPs_per_bp_plot <- ggplot(data = cluster_features_master_binned, aes(x=tag,y=SNPs_per_bp,colour=constraint_type, fill=constraint_type)) + 
  #geom_jitter(alpha=0.2) +
  geom_boxplot(alpha=0.3,outlier.shape=NA) + 
  labs(x='Intergenic cluster length (bp)', y='SNPs per bp') +
  #guides(color=FALSE) +
  scale_colour_manual(values=c("tomato","grey2")) +
  scale_fill_manual(values=c("firebrick3","grey1")) +
  theme_bw() + 
  theme(axis.text=element_text(size=10),axis.text.x=element_text(angle = 45,hjust = 1),axis.title=element_text(size=12,face="bold")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),legend.position = "none",axis.title.x = element_blank(),axis.text.x = element_blank()) +
  coord_cartesian(ylim=c(0,0.6))

SNPs_per_bp_plot

# Test for significant differences:
t.test(cluster_features_master_binned$SNPs_per_bp[constraint], cluster_features_master_binned$SNPs_per_bp[nonconstraint])
lapply(split(cluster_features_master_binned, cluster_features_master_binned$tag), function(i){plot(lm(constraint_binomial~SNPs_per_bp, data=i))})
lapply(split(cluster_features_master_binned, cluster_features_master_binned$tag), function(i){anova(lm(constraint_binomial~SNPs_per_bp, data=i))})

# Recombination rate
rr_plot <- ggplot(data = cluster_features_master_binned, aes(x=tag,y=recombination_rate,colour=constraint_type, fill=constraint_type)) + 
  #geom_jitter(alpha=0.2) +
  geom_boxplot(alpha=0.3,outlier.shape=NA) + 
  labs(x='Intergenic cluster length (bp)', y='Recombination rate (cM/Mb)') +
  #guides(color=FALSE) +
  scale_colour_manual(values=c("tomato","grey2")) +
  scale_fill_manual(values=c("firebrick3","grey1")) +
  theme_bw() + 
  theme(axis.text=element_text(size=10),axis.title=element_text(size=12,face="bold")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),legend.position = "none",axis.title.x = element_blank(),axis.text.x = element_blank()) +
  coord_cartesian(ylim=c(0,0.4))

rr_plot

# Test for significant differences:
t.test(cluster_features_master_binned$recombination_rate[constraint], cluster_features_master_binned$recombination_rate[nonconstraint])
lapply(split(cluster_features_master_binned, cluster_features_master_binned$tag), function(i){plot(lm(constraint_binomial~recombination_rate, data=i))})
lapply(split(cluster_features_master_binned, cluster_features_master_binned$tag), function(i){anova(lm(constraint_binomial~recombination_rate, data=i))})



# Combine plots
feature_plots <- grid.arrange(cluster_size_histogram, gene_distance_plot, ccre_distance_plot, dhs_distance_plot, loop_anchor_distance_plot, GWAS_SNPs_plot, av_MAF_plot, SNPs_per_bp_plot, 
                              nrow = 4, ncol=2)
# Axis aren't aligned correctly... try this:
gA <- ggplotGrob(cluster_size_histogram)
gB <- ggplotGrob(gene_distance_plot)
gC <- ggplotGrob(ccre_distance_plot)
gD <- ggplotGrob(dhs_distance_plot)
gE <- ggplotGrob(loop_anchor_distance_plot)
gF <- ggplotGrob(GWAS_SNPs_plot)
gG <- ggplotGrob(av_MAF_plot)
gH <- ggplotGrob(SNPs_per_bp_plot)
gI <- ggplotGrob(rr_plot)
maxWidth = grid::unit.pmax(gA$widths[2:5], gB$widths[2:5],gC$widths[2:5],gD$widths[2:5],gE$widths[2:5],gF$widths[2:5],gG$widths[2:5],gH$widths[2:5],gI$widths[2:5])
gA$widths[2:5] <- as.list(maxWidth)
gB$widths[2:5] <- as.list(maxWidth)
gC$widths[2:5] <- as.list(maxWidth)
gD$widths[2:5] <- as.list(maxWidth)
gE$widths[2:5] <- as.list(maxWidth)
gF$widths[2:5] <- as.list(maxWidth)
gG$widths[2:5] <- as.list(maxWidth)
gH$widths[2:5] <- as.list(maxWidth)
gI$widths[2:5] <- as.list(maxWidth)
feature_plots <- grid.arrange(gA, gB,gC,gD,gI,gE,gG,gH, nrow=4, ncol=2)


# Add labels to the arranged plots
feature_plots_labelled <- as_ggplot(feature_plots) + # transform to a ggplot
  draw_plot_label(label = c("A", "B", "C","D","E","F","G","H"), size = 15,
                  x = c(0,0.5), y = c(1,1,0.75,0.75,0.5,0.5,0.25,0.25)) # Add labels
feature_plots_labelled 


## Just plot size distribution and AFs

# Combine plots
significant_feature_plots <- grid.arrange(cluster_size_histogram, SNPs_per_bp_plot,av_MAF_plot,  
                                          nrow = 3, ncol=1)
# Axis aren't aligned correctly... try this:
gA <- ggplotGrob(cluster_size_histogram)
gB <- ggplotGrob(SNPs_per_bp_plot)
gC <- ggplotGrob(av_MAF_plot)


maxWidth = grid::unit.pmax(gA$widths[2:5],gB$widths[2:5],gC$widths[2:5])
gA$widths[2:5] <- as.list(maxWidth)
gB$widths[2:5] <- as.list(maxWidth)
gC$widths[2:5] <- as.list(maxWidth)
significant_feature_plots <- grid.arrange(gA,gB,gC, nrow=3, ncol=1)


# Add labels to the arranged plots
significant_feature_plots_labelled <- as_ggplot(significant_feature_plots) + # transform to a ggplot
  draw_plot_label(label = c("A", "B", "C"), size = 15,
                  x = c(0,0,0), y = c(1,0.66,0.33)) # Add labels
significant_feature_plots_labelled 

