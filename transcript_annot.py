#!/usr/bin/python3
import re
import argparse

gffannot = open('/uufs/chpc.utah.edu/common/home/gompert-group1/data/callosobruchus/Annotation/main_gff_ref/GCA_900659725.1_ASM90065972v1_genomic.gff','r')
transcount = open('allsamples_featurecounts.txt', 'r')
outfile = open('allsamples_featurecounts_annot.txt', 'w')

iddict = {} #region and geneid from annotation file
for line in gffannot:
	if line[0] == '#':
		continue
	line = line.split('\t')
	if line[2] == 'mRNA':
		linesub = line
		region = linesub[2]
		geneid = linesub[8]
		ids = geneid.split(';')
		tid = ids[0]
		transid1 = tid.split('=')
		transid = transid1[1]
		aid = ids[5]
		annotid1 = aid.split('=')
		annotid = annotid1[1]
		iddict[transid] = annotid

outheader = "Geneid\tChr\tStart\tEnd\tcmac16\tcmac17\tcmac18\tcmac19\tcmac28\tcmac29\tcmac31\tcmac40\tcmac41\tcmac45\tcmac58\tcmac59\t cmac63\tcmac64\tcmac65\tcmac66\tcmac67\tcmac68\tcmac69\tcmac6\tcmac70\tcmac71\tcmac72\tcmac73\tcmac7\tAnnotation\n"
outfile.write(outheader)
for line in transcount:
	line = line.strip('\n')
	line = line.split(' ')
	tids = line[0]
	w = line
	w.append('missing')
	if tids in iddict:
		ann = iddict[tids]
		w[29] = ann
		#outfile.write('\t'.join(str(w)))
		w = map(str,w)
		wline = "\t".join(w) + '\n'
		outfile.write(outheader)
		outfile.write(wline)
		
