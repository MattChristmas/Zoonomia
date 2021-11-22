library(tidyverse)
library(ggpubr)
library(gridExtra)
library(cowplot)
library(dplyr)

# This script is for calculating the proportion of constraint across the genome based on genome-wide phyloP scores
# compared to phyloP scores for ancestral repeats (nearly-neutrally evolving)

# Import cumulative frequencies of AR phylop scores

# convert to distribution table
AR_frequency_counts <- read_delim("chimp_AR_phylop_score_frequency_table_gt_-1.5.txt", 
                                  "\t", escape_double = FALSE, trim_ws = TRUE)

AR_cumulative_frequencies <- AR_frequency_counts %>% 
  arrange(phylop) %>% 
  mutate(cum_frequency=cumsum(count))
AR_cumulative_frequencies
AR_total_count=sum(AR_frequency_counts$count)
AR_cumulative_frequencies$cf = AR_cumulative_frequencies$cum_frequency/AR_total_count
AR_cumulative_frequencies$species = "chimp_AR"
write.table(AR_cumulative_frequencies, file="chimp_AR_cf_table.txt", row.names=FALSE, sep="\t",quote=FALSE)



# Use the genome-wide distribution of phylop counts over -1.5
genome_phylop_frequency_counts <- read_delim("chimp_genome_phylop_score_frequency_table_gt-1.5.txt", 
                            "\t", escape_double = FALSE, trim_ws = TRUE)

# Information on how to plot ecdf on frequency data from here: http://www.amitsharma.in/post/cumulative-distribution-plots-for-frequency-data-in-r/

cumulative_frequencies <- genome_phylop_frequency_counts %>% 
  arrange(phylop) %>% 
  mutate(cum_frequency=cumsum(count))
cumulative_frequencies

total_count=sum(cumulative_frequencies$count)
cumulative_frequencies$cf = cumulative_frequencies$cum_frequency/total_count
cumulative_frequencies$species = "chimpanzee"
write.table(cumulative_frequencies, file="chimp_cf_table.txt", row.names=FALSE, sep="\t",quote=FALSE)


# Plot ecdfs for ARs and genome together
ggplot(cumulative_frequencies, aes(x=phylop,y=cf)) +
  geom_step(colour="red") + 
  geom_step(data=AR_cumulative_frequencies, colour="grey") + 
  xlab("phyloP") + 
  ylab("Cumulative frequency") +
  theme_bw() + 
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))


# Now, calculate the value of πa, lower bound for the fraction of the genome under purifying selection

phylop_values <- cumulative_frequencies$phylop
s_list=vector(mode = "list", length = 0) #for storing the score values

for (score in phylop_values) {
  g <- subset(cumulative_frequencies, phylop==score)
  g_cf = g$cf
  ar <- subset(AR_cumulative_frequencies, phylop==score)
  
  
  if (dim(ar)[1] > 0) { # check we have a vale for this score in the ar dataframe
    ar_cf = ar$cf
    s = g_cf/ar_cf
    s_list <- c(s_list,s)  
  }

}
s_list<- unlist(s_list)
min_s <- min(s_list)

# Calculate πa, the fraction of sites in class a that are under selection

pi_a <- 1 - min_s

chimp_genome_size = 3050398082
total_phylop_scores = sum(cumulative_frequencies$count)
constraint_bases = total_phylop_scores*pi_a
prop_whole_genome = constraint_bases/chimp_genome_size

