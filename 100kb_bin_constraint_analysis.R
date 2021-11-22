library(tidyverse)
library("devtools")
#install_github("jdstorey/qvalue")
library(qvalue)
library(boot)
library(ggrepel)

options(scipen=10000)
# A function to calculate mean from the data for bootstrapping

# Analysis originally performed by Patrick F. Sullivan, script modified by Matthew J. Christmas


samplemean <- function(x, d) {  # This is the function to calculate mean from the data (x), using a bootstrap sample, d.
  return(mean(x[d]))
}

# Read in 100kb bins data frame (supplementary data file S2)
bindata_100kb <- read_delim("bindata.100000.tsv", 
                             "\t", escape_double = FALSE, trim_ws = TRUE)

b <- bindata_100kb %>%
  filter(chr != "chrY", # Remove Y chr
         zoonomiaBases>80000, # At least 80% of positions have phyloP values
         is.na(zoonomiaConsBases)==FALSE,
         (end-start)==100000) %>%    # drop 24
  mutate(zooFrac = zoonomiaConsBases/zoonomiaBases) %>%
  select(chr:end, starts_with("zoo"), k24gt90.sum, cds.distinct.sum, cCRE.sum, dhs.sum)
nrow(a)-nrow(b)

write.table(b, file="100kb_bins_analysis.txt", sep="\t",row.names = FALSE, quote=FALSE)

# Run the linear model on sqrt values
lm <- lm(sqrt(zoonomiaConsBases) ~ sqrt(k24gt90.sum) + sqrt(zoonomiaBases) + sqrt(zoonomiaNspeciesLow) + sqrt(cds.distinct.sum),data = b)
d1 <- glance(lm)   # model info, R2 etc
d2 <- tidy(lm)     # estimates, SE, P, etc
d3 <- augment(lm, data = b) %>%
  mutate(stdresidP = 2*pnorm(abs(.std.resid), lower.tail = F))

# add q-values
d4 <- qvalue(as.vector(d3$stdresidP))
e <- bind_cols(d3, tibble(d4[["qvalues"]])) %>% 
  rename(qval = `d4[[\"qvalues\"]]`) %>%
  mutate(sigBin = (qval <= 0.05))
fwrite(e, paste0(TO, "/PFS_100kb_windows.full.df.txt"), sep = "\t")

# for sigBin, combine adjacent bins
f <- e %>%
  filter(sigBin == TRUE) %>%
  mutate(.std.resid = round(.std.resid, 2),
         qval = round(qval, 6)) %>%
  select(chr, start, end, qval, .std.resid, zoonomiaConsBases, cds.distinct.sum)
fwrite(f, file = "TMPa", col.names = F, sep = "\t")
system(paste("bedtools merge -i TMPa -d 100 -c 4,5,6,7 -o min,collapse,sum,sum  > TMPb"))
g <- fread("TMPb") %>%
  mutate(len = as.integer(V3-V2),
         fracConsFdr05 = V6/len,
         fracCDS = V7/len) %>%
  select(chr=V1, start=V2, end=V3, len, qval=V4, .std.resid=V5, fracConsFdr05, fracCDS)

fwrite(g, paste0(TO, "/100kb_windows.sig.bins.txt"), sep = "\t")

# Reimport after adding gene names
sig_bins_w.genes <- read_delim("100kb_windows.sig.bins.txt", 
                           "\t", escape_double = FALSE, trim_ws = TRUE)
sig_bins_w.genes$chr_f <- as.factor(sig_bins_w.genes$chr_f)
# manhattan plot 
k <- e %>% 
  filter(sigBin==FALSE) #%>%
  mutate(chrN = case_when(chr=="chrX" ~ "23", chr=="chrY" ~ "24", TRUE ~ str_sub(chr, 4)),
         chrN = as.numeric(chrN)) %>%
  select(chrN, start, end, qval)
  #bind_rows(select(j, chrN, start, end, qval, lbl)) %>%
  #select(chrN, end, qval, lbl)

chr_labels <- as.factor(c("1","2","3","4","5","6","7","8","9","10","11",
                "12","13","14","15","16","17","18","19","20","21","22",
                "X"))
names(chr_labels) <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11",
  "chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22",
  "chrX")

chr_names = list(
  "chr1"="1",
  "chr2"="2",
  "chr3"="3",
  "chr4"="4",
  "chr5"="5",
  "chr6"="6",
  "chr7"="7",
  "chr8"="8",
  "chr9"="9",
  "chr10"="10",
  "chr11"="11",
  "chr12"="12",
  "chr13"="13",
  "chr14"="14",
  "chr15"="15",
  "chr16"="16",
  "chr17"="17",
  "chr18"="18",
  "chr19"="19",
  "chr20"="20",
  "chr21"="21",
  "chr22"="22",
  "chrX"="X"
)

chrom_labeller <- function(variable,value){
  return(chr_names[value])
}
# Make chr as factor in dfs for plotting in correct order
g$chr_f <- factor(g$chr, levels=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11",
                                  "chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22",
                                  "chrX"))
k$chr_f <- factor(k$chr, levels=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11",
                                  "chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22",
                                  "chrX"))

sig_bins_w.genes_subset <- subset(sig_bins_w.genes, qval<0.006)

# plot manhattan plot

window_constraint_man_plot <- ggplot(data=g, aes(x=start, y=-log10(qval), colour=chr)) + # g contains only significant bins
  geom_point(alpha=0.8,size=0.6) +
  geom_point(data=k, alpha=0.8,size=0.6) + # k contains only non-significant bins
  facet_grid(~chr_f, scales = "free", space="free_x", switch = "x", labeller = chrom_labeller) + # Plot chromosomes separately
  geom_hline(yintercept = -log10(0.05), colour="#30110D", size=0.5, alpha=0.4) +
  xlab("Chromosome") + 
  ylab(bquote(-log[10]("q-value"))) +
  scale_y_continuous(expand=c(0,0)) +
  #scale_y_continuous(limits = c(0, 4.8), expand = c(0, 0)) +
  geom_text_repel(data=sig_bins_w.genes_subset, aes(label = gene), size=2, na.rm = TRUE, min.segment.length = unit(0, "lines")) +
  scale_colour_manual(values=c("dodgerblue2","red2","dodgerblue2","red2","dodgerblue2","red2","dodgerblue2","red2","dodgerblue2","red2","dodgerblue2","red2","dodgerblue2","red2",
                               "dodgerblue2","red2","dodgerblue2","red2","dodgerblue2","red2","dodgerblue2","red2","dodgerblue2")) +
  theme_bw() + 
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16,face="bold")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), axis.text.x = element_blank(), 
        axis.ticks.x=element_blank(), strip.background = element_blank(),strip.placement = "outside", 
        panel.spacing=unit(0.1, "lines"), legend.position="none")

window_constraint_man_plot

