library(tidyverse)
library(ggpubr)
options(scipen=10000)

# This script is for analysing constraint in human repeats

# Read in repeats with constraint
repeats <- read_delim("human_repeats_constraint_counts.txt","\t", escape_double = FALSE, trim_ws = TRUE) 

# and without constraint
repeats_nc <- read_delim("human_repeats_no_constraint.txt","\t", escape_double = FALSE, trim_ws = TRUE) 
repeats_nc$length <- abs(repeats_nc$stop-repeats_nc$start)
repeats_nc$constraint_positions <- 0
repeats_nc$non_constraint_positions <- repeats_nc$length

repeats$prop_constraint <- repeats$constraint_positions/repeats$length
repeats_nc$prop_constraint <- 0

# combine all repeats, those with constraint sites and those without
all_repeats <- rbind(repeats,repeats_nc)
# make subsets for DNA transposons and simle repeats
DNA <- subset(all_repeats, repeat_type=="DNA")
simple_repeats <- subset(all_repeats, repeat_type=="Simple_repeat")



repeats_prop_constraint_plot <- ggplot(data=repeats, aes(x=repeat_type, y=prop_constraint)) +
  #geom_line() +
  geom_boxplot() +
  #scale_x_continuous(breaks = round(seq(min(bsyl_ID_no_tail_no_sisters_mlma$bp/1000000), max(bsyl_ID_no_tail_no_sisters_mlma$bp/1000000), by = 5),1)) +
  theme_bw() + 
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16,face="bold")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))
repeats_prop_constraint_plot


total_repeat_constraint_positions <- sum(repeats$constraint_positions)
total_repeat_non_constraint_positions <- sum(repeats$non_constraint_positions) + sum(repeats_nc$length)

# Set up df for storing values for each repeat type
repeat_constraint_df <- data.frame(type=character(),
                 total_constraint_positions=double(), 
                 constraint_prop=double(), 
                 total_non_constraint_positions=double(), 
                 non_constraint_prop=double(), 
                 diff = double(),
                 stringsAsFactors=FALSE) 


repeat_types <- unique(repeats$repeat_type)
repeat_types_nc <- unique(repeats_nc$repeat_type)
types <- c("Simple_repeat", "Low_complexity", "rRNA", "tRNA","DNA","Satellite","LINE","RC","srpRNA",
           "SINE","Retroposon","snRNA","scRNA","LTR","RNA","Unknown")

for (i in types) { # loop through all the types and get the data we need
  
  temp <- repeats$repeat_type == i
  total_constraint <- sum(repeats$constraint_positions[temp])
  constraint_prop <- total_constraint/total_repeat_constraint_positions
  
  temp_nc <- repeats_nc$repeat_type == i
  total_non_constraint <- sum(repeats$non_constraint_positions[temp]) + sum(repeats_nc$length[temp_nc])
  non_constraint_prop <- total_non_constraint/total_repeat_non_constraint_positions
  
  diff = (constraint_prop - non_constraint_prop)*100
  
  next_row <- c(i,total_constraint,constraint_prop,total_non_constraint,non_constraint_prop,diff)
  repeat_constraint_df[nrow(repeat_constraint_df) + 1, ] <- next_row

}
repeat_constraint_df$diff <- as.double(repeat_constraint_df$diff)


constraint_excess_plot <- repeat_constraint_df %>%
  #arrange(diff) %>%
  mutate(type=fct_reorder(type,-diff)) %>%
ggplot(aes(x=type, y=diff)) +
  geom_col() +
  xlab("Repeat type") +
  ylab("Constraint excess (%)") +
  theme_bw() + 
  theme(axis.text=element_text(size=12),axis.title=element_text(size=12,face="bold")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

rm(repeats)
rm(repeats_nc)
rm(repeats_prop_constraint_plot)

## Look at GC content of repeats
repeats_gc <- read_delim("./GC_content/human_repeats_GC_averages.txt","\t", escape_double = FALSE, trim_ws = TRUE) 

repeats_gc %>%
  subset(repeat_type %in% types) %>%
  mutate(repeat_type=fct_reorder(repeat_type,-mean_prop_GC)) %>%
  ggplot(aes(x=repeat_type, y=mean_prop_GC)) +
  geom_col() +
  xlab("Repeat type") +
  ylab("Mean proportion GC") +
  scale_y_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5,0.6), expand=c(0,0), limits=c(0,0.6)) +
  theme_bw() + 
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16,face="bold")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Look at correlation between constraint and GC content for DNA transposons and simple repeats
# Bring in GC data:
DNA_gc <- read_delim("./GC_content/DNA_GC_chr_renamed.txt","\t", escape_double = FALSE, trim_ws = TRUE) 
simple_repeat_gc <- read_delim("./GC_content/Simple_repeat_GC_chr_renamed.txt","\t", escape_double = FALSE, trim_ws = TRUE) 

# combine with other DFs

DNA_all <- merge(x=DNA, y=DNA_gc[,c("chrom","start","prop_GC")], by=c("chrom","start"))
simple_repeats_all <- merge(x=simple_repeats, y=simple_repeat_gc[,c("chrom","start","prop_GC")], by=c("chrom","start"))
write.table(simple_repeats_all, file="simple_repeats_constraint_GC.bed", sep="\t",quote=FALSE,row.names = FALSE)
# plots
#DNA transposons
ggplot(data=DNA_all,aes(x=prop_GC, y=prop_constraint)) +
  geom_point(size=0.4, alpha=0.5, colour="darkgrey") +
  xlab("Proportion GC") +
  ylab("Proportion constraint") +
  theme_bw() + 
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16,face="bold")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

cor.test(DNA_all$prop_GC, DNA_all$prop_constraint)

# Pearson's product-moment correlation
# 
# data:  DNA_all$prop_GC and DNA_all$prop_constraint
# t = -9.3869, df = 484000, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.01630819 -0.01067472
# sample estimates:
#         cor 
# -0.01349156

#Simple repeats
ggplot(data=simple_repeats_all,aes(x=prop_GC, y=prop_constraint)) +
  geom_point(size=0.4, alpha=0.5, colour="darkgrey") +
  xlab("Proportion GC") +
  ylab("Proportion constraint") +
  theme_bw() + 
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16,face="bold")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

cor.test(simple_repeats_all$prop_GC, simple_repeats_all$prop_constraint)


# Pearson's product-moment correlation
# 
# data:  simple_repeats_all$prop_GC and simple_repeats_all$prop_constraint
# t = 192.14, df = 678663, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.2248741 0.2293869
# sample estimates:
#       cor 
# 0.2271317



# Look at relationship between distance to gene and constraint in simple repeats

SR_gene_dist <- read_delim("simple_repeats_constraint_GC_distance_to_genes.bed","\t", escape_double = FALSE, trim_ws = TRUE) 
SR_gene_distance_plot <- SR_gene_dist %>%
  subset(distance<5000000) %>%
ggplot(aes(x=distance, y=prop_constraint)) +
  geom_point(size=0.4, alpha=0.5, colour="darkgrey") +
  xlab("Distance to nearest gene (Mb)") +
  ylab("Proportion constraint") +
  scale_x_continuous(breaks=c(0,1000000,2000000,3000000,4000000,5000000), labels=c(0,1,2,3,4,5)) +
  theme_bw() + 
  theme(axis.text=element_text(size=12),axis.title=element_text(size=12,face="bold")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))


SR_gene_dist_subset<- subset(SR_gene_dist, distance<5000000)
cor.test(SR_gene_dist_subset$distance, SR_gene_dist_subset$prop_constraint, method="spearman")

# Pearson's product-moment correlation
# 
# data:  SR_gene_dist_subset$distance and SR_gene_dist_subset$prop_constraint
# t = -60.538, df = 696502, p-value < 0.00000000000000022
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.07468357 -0.07001120
# sample estimates:
#         cor 
# -0.07234778 

# Are the majority of simple repeats near genes?
SR_gene_dist %>%
  subset(distance<5000000) %>%
  ggplot(aes(x=distance)) +
  geom_histogram(binwidth = 500000, colour="darkgrey") +
  xlab("Distance to nearest gene (Mb)") +
  #ylab("Proportion constraint") +
  #scale_x_continuous(breaks=c(0,1000000,2000000,3000000,4000000,5000000), labels=c(0,1,2,3,4,5)) +
  theme_bw() + 
  theme(axis.text=element_text(size=12),axis.title=element_text(size=12,face="bold")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))

SR_gene_dist_10kb <- subset(SR_gene_dist, distance<10000)
376344/696589

# Plot excess constraint and SR gene distance plots together



ggarrange(constraint_excess_plot,SR_gene_distance_plot, labels = c("A","B"), ncol=1,nrow=2)



