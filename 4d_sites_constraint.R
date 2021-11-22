library(tidyverse)
library(ggrepel)
library(ggpubr)

# This script is for analysing constraint at four-fold degenerate sites and overlap with ENCODE3 transcription-factor
# binding site signals (encRegTfbsClusteredWithCells.hg38.bed; available from UCSC)

#Written by Matthew Christmas, 22/11/2021, matthew.christmas aT imbim dot uu dot se

four_ds <- readRDS("fourDS_sites.rda")
four_ds$start <- four_ds$end
four_ds_reduced <- select(four_ds, chr, start, end, coding_base, amino_acid,phyloP.am, Nspecies)

# Export out to intersect with TFBSs using bedtools
write.table(four_ds_reduced, file="all_4D_sites.bed", sep="\t", quote=FALSE, row.names = FALSE)

# Bedtool intersect:

#  bedtools intersect -c -a all_4D_sites.bed -b encRegTfbsClusteredWithCells.hg38.bed > TFBS_4D_sites_intersect.bed

# Import back in

TFBS_4D_sites_intersect <- read_delim("TFBS_4D_sites_intersect.bed", 
                                      delim = "\t", escape_double = FALSE, 
                                      trim_ws = TRUE)

sum(TFBS_4D_sites_intersect$in_TFBS=="TRUE")
sum(TFBS_4D_sites_intersect$in_TFBS=="FALSE")

# Intersect with cCREs as well, import and merge
cCRE_4D_sites_intersect <- read_delim("cCRE_4D_sites_intersect.bed", 
                                     delim = "\t", escape_double = FALSE, 
                                     trim_ws = TRUE)

# And Greg Andrews' TFBSs
GA_TFBS_4D_sites_intersect <- read_delim("GA_TFBS_4D_sites_intersect.bed", 
                                      delim = "\t", escape_double = FALSE, 
                                      trim_ws = TRUE)

merge_4D_sites_intersect <- merge(TFBS_4D_sites_intersect,cCRE_4D_sites_intersect[ , c("chr","start", "cCRE_overlaps")],by=c("chr","start"))
merge_4D_sites_intersect <- merge(merge_4D_sites_intersect,GA_TFBS_4D_sites_intersect[ , c("chr","start", "GA_TFBS_overlaps")],by=c("chr","start"))
# Remove sites without phyloP values
merge_4D_sites_intersect <- merge_4D_sites_intersect[complete.cases(merge_4D_sites_intersect[,6]),]
# Remove duplicated rows
merge_4D_sites_intersect_unique <- 
  merge_4D_sites_intersect %>%
    distinct()

merge_4D_sites_intersect_unique$constrained <- "TRUE"
non_constrained <- merge_4D_sites_intersect_unique$phyloP.am <=2.270
merge_4D_sites_intersect_unique$constrained[non_constrained] <- "FALSE"
merge_4D_sites_intersect_unique$in_TFBS <- "FALSE"
in_TFBS <- merge_4D_sites_intersect_unique$TFBS_overlaps>0
merge_4D_sites_intersect_unique$in_TFBS[in_TFBS] <- "TRUE"

merge_4D_sites_intersect_unique$in_cCRE <- "FALSE"
in_cCRE <- merge_4D_sites_intersect_unique$cCRE_overlaps>0
merge_4D_sites_intersect_unique$in_cCRE[in_cCRE] <- "TRUE"


merge_4D_sites_intersect_unique$in_GA_TFBS <- "FALSE"
in_GA_TFBS <- merge_4D_sites_intersect_unique$GA_TFBS_overlaps>0
merge_4D_sites_intersect_unique$in_GA_TFBS[in_GA_TFBS] <- "TRUE"

TFBS_sites <- subset(merge_4D_sites_intersect_unique, in_TFBS=="TRUE")
nonTFBS_sites <- subset(merge_4D_sites_intersect_unique, in_TFBS=="FALSE")
mean(TFBS_sites$phyloP.am,na.rm=TRUE)
mean(nonTFBS_sites$phyloP.am, na.rm=TRUE)
#t.test(TFBS_sites$phyloP.am, nonTFBS_sites$phyloP.am)


cCRE_sites <- subset(merge_4D_sites_intersect_unique, in_cCRE=="TRUE")
noncCRE_sites <- subset(merge_4D_sites_intersect_unique, in_cCRE=="FALSE")
mean(cCRE_sites$phyloP.am,na.rm=TRUE)
mean(noncCRE_sites$phyloP.am, na.rm=TRUE)
#t.test(cCRE_sites$phyloP.am, noncCRE_sites$phyloP.am)


GA_TFBS_sites <- subset(merge_4D_sites_intersect_unique, in_GA_TFBS=="TRUE")
nonGA_TFBS_sites <- subset(merge_4D_sites_intersect_unique, in_GA_TFBS=="FALSE")
mean(GA_TFBS_sites$phyloP.am,na.rm=TRUE)
mean(nonGA_TFBS_sites$phyloP.am, na.rm=TRUE)
#t.test(GA_TFBS_sites$phyloP.am, nonGA_TFBS_sites$phyloP.am)

in_annotation_sites <- subset(merge_4D_sites_intersect_unique, in_GA_TFBS=="TRUE" | in_TFBS=="TRUE" | in_cCRE=="TRUE")
not_in_annotation_sites <- subset(merge_4D_sites_intersect_unique, in_GA_TFBS=="FALSE" & in_TFBS=="FALSE" & in_cCRE=="FALSE")
mean(in_annotation_sites$phyloP.am,na.rm=TRUE)
mean(not_in_annotation_sites$phyloP.am, na.rm=TRUE)
t.test(in_annotation_sites$phyloP.am, not_in_annotation_sites$phyloP.am)

colors <- c("In TFBS" = "tomato", "Outside TFBS" = "blue")

ggplot(data=in_annotation_sites) +
  geom_histogram(aes(y=..count../sum(..count..), x=phyloP.am, fill="In TFBS"), binwidth = 1, alpha=0.5) +
  geom_histogram(data=not_in_annotation_sites, aes(y=..count../sum(..count..), x=phyloP.am, fill="Outside TFBS"), binwidth = 1, alpha=0.3) +
  labs(y="Proportion of 4D sites", x="PhyloP", fill="Location") +
  scale_fill_manual(values = colors) +
  geom_vline(xintercept = 2.270, lty=2, colour="red") +
  theme_bw() + 
  theme(axis.text=element_text(size=12),axis.title=element_text(size=16,face="bold")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), legend.position=c(0.2,0.8))



# Test for enrichment of constrained 4D sites within annotations
total_sites <- nrow(merge_4D_sites_intersect_unique)
total_constrained_sites <- sum(merge_4D_sites_intersect_unique$constrained=="TRUE")
total_nonconstrained_sites <- sum(merge_4D_sites_intersect_unique$constrained=="FALSE")
total_sites_in_annotations <- nrow(in_annotation_sites)
total_constrained_sites_in_annotations <- sum(in_annotation_sites$constrained=="TRUE")
total_nonconstrained_sites_in_annotations <- sum(in_annotation_sites$constrained=="FALSE")

contingency_table <- data.frame("total"=c(total_sites,total_constrained_sites), 
                                "in_annotation"=c(total_sites_in_annotations,total_constrained_sites_in_annotations))
row.names(contingency_table)=c("all","constraint")

fisher.test(contingency_table)
chisq.test(contingency_table)



perc_all_sites_constrained <- 965406/5025952*100 # 19.2%
perc_annotated_sites_constrained <- 625728/2980106*100 # 21%

perc_all_sites_annotated <- total_sites_in_annotations/total_sites*100 # 59.29%
perc_constrained_sites_annotated <- total_constrained_sites_in_annotations/total_constrained_sites*100 # 64%
perc_nonconstrained_sites_annotated <- total_nonconstrained_sites_in_annotations/total_nonconstrained_sites*100 # 57.98%


### Look at whether certain TFs contain a surplus of constrained 4D sites
# This is all TFBSs, regardless of cluster score

TF_data <- read_delim("TF_4D_sites_constraint_counts.txt", 
                      delim = "\t", escape_double = FALSE, 
                      trim_ws = TRUE)

ggplot(data=TF_data, aes(x=prop_constraint_sites)) +
  geom_histogram(binwidth=0.01) +
  xlab("Proportion of 4D sites that are constrained") +
  theme_bw() + 
  theme(axis.text=element_text(size=12),axis.title=element_text(size=16,face="bold")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))


# use score threshold of 1000

TF_data_1000 <- read_delim("TF_4D_sites_constraint_counts_score_threshold_1000.txt", 
                      delim = "\t", escape_double = FALSE, 
                      trim_ws = TRUE)

total_constrained_4DS_overlapping_TFBS_1000 <- 209958
total_non_constrained_4DS_overlapping_TFBS_1000 <- 685499

TF_data_1000$prop_constrained <- (TF_data_1000$constraint_sites/total_constrained_4DS_overlapping_TFBS_1000)
TF_data_1000$prop_non_constrained <- (TF_data_1000$non_constraint_sites/total_non_constrained_4DS_overlapping_TFBS_1000)
TF_data_1000$perc_difference <- (TF_data_1000$prop_constrained-TF_data_1000$prop_non_constrained)*100
  
ggplot(data=TF_data_1000, aes(x=perc_difference)) +
  geom_histogram(binwidth=0.05) +
  xlab("Constrained versus non-constrained 4D sites") +
  theme_bw() + 
  theme(axis.text=element_text(size=12),axis.title=element_text(size=16,face="bold")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))

# Define outliers for plotting
is_outlier <- function(x) {
  return(x >= quantile(x, 0.97) | x <= quantile(x, 0.03))
}


# Create a forest plot, with TFs ordered in descending order of perc_difference
plotA <-TF_data_1000 %>%
  mutate(is_outlier=ifelse(is_outlier(perc_difference), TF, NA)) %>%
  #mutate(is_CTCF=ifelse(TF=="CTCF", TF, NA)) %>%
  arrange(desc(perc_difference)) %>%
  mutate(TF=fct_reorder(TF,perc_difference)) %>%
ggplot(aes(x=perc_difference, y=TF)) +
  geom_point(size=0.3) +
  geom_vline(xintercept = -0.3847068, lty=2, colour="darkgrey") + # add lines to indicate 95th percentile 
  geom_vline(xintercept =0.3518533, lty=2, colour="darkgrey") +
  #geom_text(aes(label=is_outlier),na.rm=TRUE,position=position_jitter(width=0.2,height=0),size=2) +
  geom_text_repel(aes(label=is_outlier),na.rm=TRUE,size=3) +
  xlab("Excess constrained 4D sites in TFBS (%)") +
  scale_y_discrete(expand=c(0.01,0.01)) +
  scale_x_continuous(breaks=c(-1,0,1,2,3,4)) +
  theme_bw() + 
  theme(axis.text=element_text(size=10),axis.title=element_text(size=16,face="bold")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), axis.text.y = element_blank(),axis.ticks.y = element_blank())
plotA

quantile(TF_data_1000$perc_difference, c(.95,.97,.99))
TFs_with_constrained_TFBS <- subset(TF_data_1000, perc_difference>0.805)

write.table(TF_data_1000, file="TF_4DS_constraint_threshold_1000_proportions.txt", sep="\t",quote=FALSE, row.names = FALSE)

# use score threshold of 800

TF_data_800 <- read_delim("TF_4D_sites_constraint_counts_score_threshold_800.txt", 
                           delim = "\t", escape_double = FALSE, 
                           trim_ws = TRUE)

total_constrained_4DS_overlapping_TFBS_800 <- 243722
total_non_constrained_4DS_overlapping_TFBS_800 <- 812381

TF_data_800$prop_constrained <- (TF_data_800$constraint_sites/total_constrained_4DS_overlapping_TFBS_800)
TF_data_800$prop_non_constrained <- (TF_data_800$non_constraint_sites/total_non_constrained_4DS_overlapping_TFBS_800)
TF_data_800$perc_difference <- (TF_data_800$prop_constrained-TF_data_800$prop_non_constrained)*100

ggplot(data=TF_data_800, aes(x=perc_difference)) +
  geom_histogram(binwidth=0.05) +
  xlab("Constrained versus non-constrained 4D sites") +
  theme_bw() + 
  theme(axis.text=element_text(size=12),axis.title=element_text(size=16,face="bold")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))

write.table(TF_data_800, file="TF_4DS_constraint_threshold_800_proportions.txt", sep="\t",quote=FALSE, row.names = FALSE)
quantile(TF_data_800$prop_constraint_sites, c(.90,.95,.99))


# Look at overlap of TFBS and exons - do TFs with high-constraint TFBSs bind
# mostly to coding sequence?

TFBS_exon_intersect <- read_delim("TFBS_exon_overlap_counts.txt", 
                                      delim = "\t", escape_double = FALSE, 
                                      trim_ws = TRUE)

TFBS_exon_intersect$prop_exon <- TFBS_exon_intersect$TFBS_overlapping_exons/(TFBS_exon_intersect$TFBS_overlapping_exons+TFBS_exon_intersect$TFBS_not_overlapping_exons)


# Look at relationship between overlap with exons and 4D site constraint
TFBS_exon_intersect_plus_constraint <- merge(TFBS_exon_intersect,TF_data_1000[, c("TF","perc_difference")], by = "TF")

ggplot(data=TFBS_exon_intersect_plus_constraint, aes(x=perc_difference, y=prop_exon)) +
  geom_point() +
  #xlab("") +
  theme_bw() + 
  theme(axis.text=element_text(size=12),axis.title=element_text(size=16,face="bold")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))

# Look at only TFBS with score = 1000 (highest confidence)

TFBS_1000_exon_intersect <- read_delim("TFBS_1000_exon_overlap_counts.txt", 
                                  delim = "\t", escape_double = FALSE, 
                                  trim_ws = TRUE)

TFBS_1000_exon_intersect$prop_exon <- TFBS_1000_exon_intersect$TFBS_overlapping_exons/(TFBS_1000_exon_intersect$TFBS_overlapping_exons+TFBS_1000_exon_intersect$TFBS_not_overlapping_exons)

# Look at relationship between overlap with exons and 4D site constraint
TFBS_1000_exon_intersect_plus_constraint <- merge(TFBS_1000_exon_intersect,TF_data_1000[, c("TF","perc_difference")], by = "TF")
plotB <-TFBS_1000_exon_intersect_plus_constraint %>%
  mutate(is_outlier=ifelse(is_outlier(perc_difference), TF, NA)) %>%
  mutate(outliers=ifelse(is_outlier(perc_difference), "Outlier", NA)) %>% # for colouring outlier points
  arrange(desc(perc_difference)) %>%
  #mutate(TF=fct_reorder(TF,perc_difference)) %>%
ggplot(aes(y=perc_difference, x=prop_exon)) +
  geom_point(aes(colour=outliers), size=0.5) +
  geom_text_repel(aes(label=is_outlier),na.rm=TRUE,size=4) +
  ylab("Excess constrained 4D sites in TFBS (%)") +
  xlab("Proportion of TFBSs in exons") +
  #scale_colour_manual(values=c("darkgrey","red")) +
  scale_y_continuous(breaks=c(-1,0,1,2,3,4)) +
  theme_bw() + 
  theme(axis.text=element_text(size=12),axis.title=element_text(size=16,face="bold")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), legend.position = "none")


ggarrange(plotA, plotB)
