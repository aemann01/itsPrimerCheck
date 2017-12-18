#!/usr/bin/python3

'''Useage: python itsPrimerCheck.py input_R1.fastq input_R2.fastq sampleName'''

import sys
import gzip
from Bio import SeqIO
from itertools import zip_longest

primer_its1f = "TCCGTAGGTGAACCTGCGG"
primer_its2r = "GCTGCGTTCTTCATCGATGC"
primer_its3f = "GCATCGATGAAGAACGCAGC"
primer_its4r = "TCCTCCGCTTATTGATATGC"

r1_its1 = []
r1_its2 = []
r2_its1 = []
r2_its2 = []

readiter_R1 = SeqIO.parse(open(sys.argv[1]), "fastq")
readiter_R2 = SeqIO.parse(open(sys.argv[2]), "fastq")

for rec1, rec2 in zip_longest(readiter_R1, readiter_R2):
	#check for normal configuration
	if primer_its1f in rec1.seq[0:24] and primer_its2r in rec2.seq[0:24]:
		r1_its1.append(rec1)
		r2_its1.append(rec2)
	if primer_its3f in rec1.seq[0:24] and primer_its4r in rec2.seq[0:24]:
		r1_its2.append(rec1)
		r2_its2.append(rec2)
	#check for opposite configuration
	if primer_its2r in rec1.seq[0:24] and primer_its1f in rec2.seq[0:24]:
		r1_its1.append(rec2)
		r2_its1.append(rec1)
	if primer_its4r in rec1.seq[0:24] and primer_its3f in rec2.seq[0:24]:
		r1_its2.append(rec2)
		r2_its2.append(rec1)
	#check for one of the seqs not having primer seq, keep both
	if primer_its1f in rec1.seq[0:24] and primer_its2r not in rec2.seq[0:24]:
		r1_its1.append(rec1)
		r2_its1.append(rec2)
	if primer_its1f in rec2.seq[0:24] and primer_its2r not in rec1.seq[0:24]:
		r1_its1.append(rec2)
		r2_its1.append(rec1)
	if primer_its3f in rec1.seq[0:24] and primer_its4r not in rec2.seq[0:24]:
		r1_its2.append(rec1)
		r2_its2.append(rec2)
	if primer_its3f in rec2.seq[0:24] and primer_its4r not in rec1.seq[0:24]:
		r1_its2.append(rec2)
		r2_its2.append(rec1)

with open('%s_ITS1_R1.fastq' % sys.argv[3], "w") as outITS1_R1:
	SeqIO.write(r1_its1, outITS1_R1, "fastq")
outITS1_R1.close()

with open("%s_ITS1_R2.fastq" % sys.argv[3], "w") as outITS1_R2:
	SeqIO.write(r2_its1, outITS1_R2, "fastq")
outITS1_R2.close()

with open("%s_ITS2_R1.fastq" % sys.argv[3], "w") as outITS2_R1:
	SeqIO.write(r1_its2, outITS2_R1, "fastq")
outITS2_R1.close()

with open("%s_ITS2_R2.fastq" % sys.argv[3], "w") as outITS2_R2:
	SeqIO.write(r2_its2, outITS2_R2, "fastq")
outITS2_R2.close()


