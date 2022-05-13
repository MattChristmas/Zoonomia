#! /usr/bin/env python
#usage: script.py [fasta of query protein] [exonerate output]
from Bio import SeqIO
import sys
def main(argv):

	protlen = 0
	for seq in SeqIO.parse(sys.argv[1], 'fasta'):
		protlen = len(str(seq.seq))
	
	pairs = {}
	for line in open(sys.argv[2]):
		line = line.strip('\n')
		qstart = int(line.split('\t')[3])
		qend = line.split('\t')[4]
		if qstart in pairs:
			pairs[qstart].append(line)
		else:
			pairs[qstart] = [line]
	gaps = []

	sorted_pairs = sorted(pairs.keys())
	start = 1
	end = 1
	for s in sorted_pairs:
		for p in pairs[s]:
			qstart = int(p.split('\t')[3])
			qend = int(p.split('\t')[4])
			
			if qstart <= (end + 1):
				#no gap
				start = qstart
				if qend > end:
					end = qend
			else:
				#gap
				if start == 1:
					gaps.append([end, qstart, "start"])
				else:
					gaps.append([end, qstart, "middle"])
				start = qstart
				if qend > end:
					end = qend
	
	if end < protlen:
		gaps.append([end, protlen, 'end'])
	
	for g in gaps:
		print("gap\t" + str(g[0]) + '\t' + str(g[1]) + '\t' + str(g[1] - g[0]) + '\t' + g[2])
		


if __name__ == "__main__":
  main(sys.argv)
