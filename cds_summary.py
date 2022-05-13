#! /usr/bin/env python
#workflow: this script is run after bedtools coverage is run for all alignment files in a directory
#all bedtools coverage output files are assumed to have the suffix "_CMAH_CDS_coverage.txt"
from __future__ import division
import sys,os
def main(argv):

	print("Source\tRef\tIntact_CDS\tMissing_CDS\tTotal_CDS\tPercent_missing_CDS\tWhich_CDS_missing\tMissing_CDS_coverage")
	for f in os.listdir(os.getcwd()):
		if "_CMAH_CDS_coverage.txt" in f:
			#print(f)
			exp = f.replace("_CMAH_CDS_coverage.txt", "")
			uppers = [i for i in range(0, len(exp)) if exp[i].isupper()]
			second = uppers[1]
			source = exp[0:second].strip('_')
			ref = exp[second:]
	
			#source = "_".join(exp.split('_')[0:2])
			#ref = "_".join(exp.split('_')[2:])
			cdses = []
			missing = {}
			total_loss_bases = 0
			intact = []
			count = 0
			for line in open(f):
				line = line.strip('\n')
				cds = line.split('\t')[0] + ':' + line.split('\t')[1] + '-' + line.split('\t')[2]
				cdses.append(cds)
				coverage = line.split('\t')[-1]
				count += 1
				if float(coverage) < 1:
					missing[str(count)] = coverage
					total = int(line.split('\t')[5]) - int(line.split('\t')[4])
					total_loss_bases += total
				else:
					intact.append(str(count))

			missing_count = len(list(missing.keys()))
			intact_count = len(intact)
			total = len(cdses)
			missing_exons = sorted(list(missing.keys()))
			missing_coverages = [missing[c] for c in missing_exons]
			percent_missing = missing_count / total
			print(source + '\t' + ref + '\t' + str(intact_count) + '\t' + str(missing_count) + '\t' + str(total) + '\t' + str(percent_missing) + '\t' + ",".join(missing_exons) + '\t' + ','.join(missing_coverages) + '\t' + str(total_loss_bases))
			

if __name__ == "__main__":
  main(sys.argv)
