# Zoonomia
Scripts used in the analysis of the 241-mammal Zoonomia dataset for the upcoming publication: "Evolutionary constraint and innovation across hundreds of placental mammals" UNDER REVIEW

These scripts and data files are provided with the main purpose of reproducibility of the analyses presented in the paper. Any questions regarding the content and use of this repository should be directed to matthew dot christmas aT imbim dot uu dot se

# Estimating genome-wide constraint:
https://github.com/MattChristmas/Zoonomia/blob/3e6e0665e5e7deb9e13c1c5848f50cb425037b26/calculate_genome_constraint_fraction.R

This script provides an example of how we calculated a lower-bound estimate of the proportion of the genome under purifying selection based on phyloP values, for human, chimpanzee, house mouse, dog and bat

# Constraint at four-fold degenerate sites:
https://github.com/MattChristmas/Zoonomia/blob/a7353d2271d4e346f83ccd2c3369783d131d267b/4d_sites_constraint.R

This script was used to analyse constraint at 4-fold degenerate sites and their overlap with transcription factor binding sites

# 100kb bins constraint analysis
https://github.com/MattChristmas/Zoonomia/blob/02b026d082ebbe03a59388c6b8b38bcadb90a672/100kb_bin_constraint_analysis.R

This script was used to analyse supplementary data file 2, for measuring constraint across the human genome in 100 kb bins

# Gene desert analysis
https://github.com/MattChristmas/Zoonomia/blob/02b026d082ebbe03a59388c6b8b38bcadb90a672/Gene_deserts_constraint_analysis.R

This script was used for analysing constraint in gene deserts and relating this to the location of gene deserts near developmental transcription factors

# UNICORNs
https://github.com/MattChristmas/Zoonomia/blob/ce35c2eee7eda1a2dff15580abed8741d0fc6e44/UNICORN_analysis.R

This script was used for analysing features of UNICORNs, in particular for showing they contain less variation with lower allele frequencies than other unannotated intergenic regions

https://github.com/MattChristmas/Zoonomia/blob/main/getUnicornOverlapsWithBrainDatasets.sh

This script was used for evaluating the overlap of UNICORNs with open chromatin regions from different brain regions, different motor cortex cell types, and different brain developmental stages

# Constraint in repeats
https://github.com/MattChristmas/Zoonomia/blob/ce35c2eee7eda1a2dff15580abed8741d0fc6e44/repeats_constraint_analysis.R

This script was used for analysing the distribution of constraint across human repeat families

# Olfaction evolution
https://github.com/MattChristmas/Zoonomia/blob/main/extract_OR_hits.pl

This script was used to extract regions of the genome containing putative olfactory receptor sequences identified using tblastx
