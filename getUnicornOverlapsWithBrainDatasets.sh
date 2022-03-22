# Get genome coverage
sort -u -k1,1 -k2,2n -k3,3n /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/UNICORNs/UNICORNs.bed | genomeCoverageBed -i stdin -g /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/HumanGenome/hg38AndPatch11.chrom.sizes > /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/UNICORNs/UNICORNs_hg38Coverage.txt
sort -u -k1,1 -k2,2n -k3,3n /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/UNICORNs/UNICORNs.bed | awk 'BEGIN{OFS="\t"} $3-$2 >= 500' | genomeCoverageBed -i stdin -g /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/HumanGenome/hg38AndPatch11.chrom.sizes > /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/UNICORNs/UNICORNs_atLeast500bp_hg38Coverage.txt

# Get amount covered by Fullard et al. neuron peaks
coverageBed -a /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/UNICORNs/UNICORNs.bed -b /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/HumanDNase/ACC/cromwell-executions/atac/ec111d9d-0ad5-46b3-a9f6-edf2e15cf560/call-reproducibility_idr/execution/optimal_peak.narrowPeak.gz /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/HumanDNase/DLPFC/cromwell-executions/atac/056f86bc-e9fd-4651-972e-86378b0dc4c1/call-reproducibility_idr/execution/optimal_peak.narrowPeak.gz /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/HumanDNase/INS/cromwell-executions/atac/f772b81d-87b3-4fed-bf05-86fff10cc7d7/call-reproducibility_idr/execution/optimal_peak.narrowPeak.gz /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/HumanDNase/ITC/cromwell-executions/atac/ac9629fb-ea84-42b1-a9df-48d37ec6e4a2/call-reproducibility_idr/execution/optimal_peak.narrowPeak.gz /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/HumanDNase/MDT/cromwell-executions/atac/8e947865-c610-44cf-9ca5-bf04db58d45f/call-reproducibility_idr/execution/optimal_peak.narrowPeak.gz /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/HumanDNase/OFC/cromwell-executions/atac/d57e29ef-dfa3-47da-9f19-cfb752b72338/call-reproducibility_idr/execution/optimal_peak.narrowPeak.gz /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/HumanDNase/PMC/cromwell-executions/atac/45c4bd8c-c055-4cea-80c9-5ec0eeaf8024/call-reproducibility_idr/execution/optimal_peak.narrowPeak.gz /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/HumanDNase/PUT/cromwell-executions/atac/cee5204e-3e9a-4437-868e-9fc8f05312af/call-reproducibility_idr/execution/optimal_peak.narrowPeak.gz /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/HumanDNase/PVC/cromwell-executions/atac/15683be5-6c10-47f7-a394-2e89f16a0b27/call-reproducibility_idr/execution/optimal_peak.narrowPeak.gz /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/HumanDNase/STC/cromwell-executions/atac/60aa17f1-e5f0-49ab-926f-26f54aecf428/call-reproducibility_idr/execution/optimal_peak.narrowPeak.gz /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/HumanDNase/VLPFC/cromwell-executions/atac/452c953a-2d83-4175-8f99-3941986d4168/call-reproducibility_idr/execution/optimal_peak.narrowPeak.gz /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/HumanDNase/AMY/cromwell-executions/atac/782bcba0-574b-4818-a393-2b4d70adfe15/call-reproducibility_idr/execution/optimal_peak.narrowPeak.gz /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/HumanDNase/HIPP/cromwell-executions/atac/62febcc1-d4d9-4453-ac85-23f18541075a/call-reproducibility_idr/execution/optimal_peak.narrowPeak.gz /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/HumanDNase/NAC/cromwell-executions/atac/f47a76c2-4f55-458f-a1ac-08713d0f0608/call-reproducibility_idr/execution/optimal_peak.narrowPeak.gz > /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/UNICORNs/UNICORNs_FullardNeuNCoverage.txt

# Get amount covered by Bakken et al. neuron peaks
coverageBed -a /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/UNICORNs/UNICORNs.bed -b /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/HumanDNase/M1SnareSeq/ArchR_BICCN_huMOp_SNARE-Seq2_moreCells/BICCN_huMOp_SNARE-Seq2_moreCells_L2.3.IT_reproduciblePeaks.bed /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/HumanDNase/M1SnareSeq/ArchR_BICCN_huMOp_SNARE-Seq2_moreCells/BICCN_huMOp_SNARE-Seq2_moreCells_L5.6.NP_reproduciblePeaks.bed /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/HumanDNase/M1SnareSeq/ArchR_BICCN_huMOp_SNARE-Seq2_moreCells/BICCN_huMOp_SNARE-Seq2_moreCells_L5.ET_reproduciblePeaks.bed /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/HumanDNase/M1SnareSeq/ArchR_BICCN_huMOp_SNARE-Seq2_moreCells/BICCN_huMOp_SNARE-Seq2_moreCells_L5.IT_reproduciblePeaks.bed /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/HumanDNase/M1SnareSeq/ArchR_BICCN_huMOp_SNARE-Seq2_moreCells/BICCN_huMOp_SNARE-Seq2_moreCells_L6b_reproduciblePeaks.bed /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/HumanDNase/M1SnareSeq/ArchR_BICCN_huMOp_SNARE-Seq2_moreCells/BICCN_huMOp_SNARE-Seq2_moreCells_L6.CT_reproduciblePeaks.bed /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/HumanDNase/M1SnareSeq/ArchR_BICCN_huMOp_SNARE-Seq2_moreCells/BICCN_huMOp_SNARE-Seq2_moreCells_L6.IT.Car3_reproduciblePeaks.bed /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/HumanDNase/M1SnareSeq/ArchR_BICCN_huMOp_SNARE-Seq2_moreCells/BICCN_huMOp_SNARE-Seq2_moreCells_L6.IT_reproduciblePeaks.bed /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/HumanDNase/M1SnareSeq/ArchR_BICCN_huMOp_SNARE-Seq2_moreCells/BICCN_huMOp_SNARE-Seq2_moreCells_LAMP5_reproduciblePeaks.bed /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/HumanDNase/M1SnareSeq/ArchR_BICCN_huMOp_SNARE-Seq2_moreCells/BICCN_huMOp_SNARE-Seq2_moreCells_PVALB_reproduciblePeaks.bed /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/HumanDNase/M1SnareSeq/ArchR_BICCN_huMOp_SNARE-Seq2_moreCells/BICCN_huMOp_SNARE-Seq2_moreCells_SNCG_reproduciblePeaks.bed /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/HumanDNase/M1SnareSeq/ArchR_BICCN_huMOp_SNARE-Seq2_moreCells/BICCN_huMOp_SNARE-Seq2_moreCells_SST.CHODL_reproduciblePeaks.bed /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/HumanDNase/M1SnareSeq/ArchR_BICCN_huMOp_SNARE-Seq2_moreCells/BICCN_huMOp_SNARE-Seq2_moreCells_SST_reproduciblePeaks.bed /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/HumanDNase/M1SnareSeq/ArchR_BICCN_huMOp_SNARE-Seq2_moreCells/BICCN_huMOp_SNARE-Seq2_moreCells_VIP_reproduciblePeaks.bed > /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/UNICORNs/UNICORNs_BakkenNeuronCoverage.txt

# Get amount covered by Markenscoff-Papadimitriou et al. brain development peaks
coverageBed -a /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/UNICORNs/UNICORNs.bed -b /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/HumanDNase/TelencephalonDevelopment/GSE149268_annotation-ocr-hg38.bed.gz > /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/UNICORNs/UNICORNs_TelencephalonDevelopmentCoverage.txt

# Get amount covered by all 3 brain peak datasets
coverageBed -a /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/UNICORNs/UNICORNs.bed -b /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/HumanDNase/ACC/cromwell-executions/atac/ec111d9d-0ad5-46b3-a9f6-edf2e15cf560/call-reproducibility_idr/execution/optimal_peak.narrowPeak.gz /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/HumanDNase/DLPFC/cromwell-executions/atac/056f86bc-e9fd-4651-972e-86378b0dc4c1/call-reproducibility_idr/execution/optimal_peak.narrowPeak.gz /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/HumanDNase/INS/cromwell-executions/atac/f772b81d-87b3-4fed-bf05-86fff10cc7d7/call-reproducibility_idr/execution/optimal_peak.narrowPeak.gz /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/HumanDNase/ITC/cromwell-executions/atac/ac9629fb-ea84-42b1-a9df-48d37ec6e4a2/call-reproducibility_idr/execution/optimal_peak.narrowPeak.gz /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/HumanDNase/MDT/cromwell-executions/atac/8e947865-c610-44cf-9ca5-bf04db58d45f/call-reproducibility_idr/execution/optimal_peak.narrowPeak.gz /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/HumanDNase/OFC/cromwell-executions/atac/d57e29ef-dfa3-47da-9f19-cfb752b72338/call-reproducibility_idr/execution/optimal_peak.narrowPeak.gz /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/HumanDNase/PMC/cromwell-executions/atac/45c4bd8c-c055-4cea-80c9-5ec0eeaf8024/call-reproducibility_idr/execution/optimal_peak.narrowPeak.gz /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/HumanDNase/PUT/cromwell-executions/atac/cee5204e-3e9a-4437-868e-9fc8f05312af/call-reproducibility_idr/execution/optimal_peak.narrowPeak.gz /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/HumanDNase/PVC/cromwell-executions/atac/15683be5-6c10-47f7-a394-2e89f16a0b27/call-reproducibility_idr/execution/optimal_peak.narrowPeak.gz /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/HumanDNase/STC/cromwell-executions/atac/60aa17f1-e5f0-49ab-926f-26f54aecf428/call-reproducibility_idr/execution/optimal_peak.narrowPeak.gz /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/HumanDNase/VLPFC/cromwell-executions/atac/452c953a-2d83-4175-8f99-3941986d4168/call-reproducibility_idr/execution/optimal_peak.narrowPeak.gz /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/HumanDNase/AMY/cromwell-executions/atac/782bcba0-574b-4818-a393-2b4d70adfe15/call-reproducibility_idr/execution/optimal_peak.narrowPeak.gz /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/HumanDNase/HIPP/cromwell-executions/atac/62febcc1-d4d9-4453-ac85-23f18541075a/call-reproducibility_idr/execution/optimal_peak.narrowPeak.gz /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/HumanDNase/NAC/cromwell-executions/atac/f47a76c2-4f55-458f-a1ac-08713d0f0608/call-reproducibility_idr/execution/optimal_peak.narrowPeak.gz /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/HumanDNase/M1SnareSeq/ArchR_BICCN_huMOp_SNARE-Seq2_moreCells/BICCN_huMOp_SNARE-Seq2_moreCells_L2.3.IT_reproduciblePeaks.bed /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/HumanDNase/M1SnareSeq/ArchR_BICCN_huMOp_SNARE-Seq2_moreCells/BICCN_huMOp_SNARE-Seq2_moreCells_L5.6.NP_reproduciblePeaks.bed /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/HumanDNase/M1SnareSeq/ArchR_BICCN_huMOp_SNARE-Seq2_moreCells/BICCN_huMOp_SNARE-Seq2_moreCells_L5.ET_reproduciblePeaks.bed /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/HumanDNase/M1SnareSeq/ArchR_BICCN_huMOp_SNARE-Seq2_moreCells/BICCN_huMOp_SNARE-Seq2_moreCells_L5.IT_reproduciblePeaks.bed /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/HumanDNase/M1SnareSeq/ArchR_BICCN_huMOp_SNARE-Seq2_moreCells/BICCN_huMOp_SNARE-Seq2_moreCells_L6b_reproduciblePeaks.bed /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/HumanDNase/M1SnareSeq/ArchR_BICCN_huMOp_SNARE-Seq2_moreCells/BICCN_huMOp_SNARE-Seq2_moreCells_L6.CT_reproduciblePeaks.bed /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/HumanDNase/M1SnareSeq/ArchR_BICCN_huMOp_SNARE-Seq2_moreCells/BICCN_huMOp_SNARE-Seq2_moreCells_L6.IT.Car3_reproduciblePeaks.bed /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/HumanDNase/M1SnareSeq/ArchR_BICCN_huMOp_SNARE-Seq2_moreCells/BICCN_huMOp_SNARE-Seq2_moreCells_L6.IT_reproduciblePeaks.bed /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/HumanDNase/M1SnareSeq/ArchR_BICCN_huMOp_SNARE-Seq2_moreCells/BICCN_huMOp_SNARE-Seq2_moreCells_LAMP5_reproduciblePeaks.bed /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/HumanDNase/M1SnareSeq/ArchR_BICCN_huMOp_SNARE-Seq2_moreCells/BICCN_huMOp_SNARE-Seq2_moreCells_PVALB_reproduciblePeaks.bed /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/HumanDNase/M1SnareSeq/ArchR_BICCN_huMOp_SNARE-Seq2_moreCells/BICCN_huMOp_SNARE-Seq2_moreCells_SNCG_reproduciblePeaks.bed /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/HumanDNase/M1SnareSeq/ArchR_BICCN_huMOp_SNARE-Seq2_moreCells/BICCN_huMOp_SNARE-Seq2_moreCells_SST.CHODL_reproduciblePeaks.bed /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/HumanDNase/M1SnareSeq/ArchR_BICCN_huMOp_SNARE-Seq2_moreCells/BICCN_huMOp_SNARE-Seq2_moreCells_SST_reproduciblePeaks.bed /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/HumanDNase/M1SnareSeq/ArchR_BICCN_huMOp_SNARE-Seq2_moreCells/BICCN_huMOp_SNARE-Seq2_moreCells_VIP_reproduciblePeaks.bed /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/HumanDNase/TelencephalonDevelopment/GSE149268_annotation-ocr-hg38.bed.gz > /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/UNICORNs/UNICORNs_BrainRegionCellTypeDevelopmentCoverage.txt