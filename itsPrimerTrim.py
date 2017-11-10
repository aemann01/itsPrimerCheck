#!/usr/bin/python3

from Bio import SeqIO
from itertools import izip_longest

primer_its1f = "TCCGTAGGTGAACCTGCGG"
primer_its2r = "GCTGCGTTCTTCATCGATGC"
primer_its3f = "GCATCGATGAAGAACGCAGC"
primer_its4r = "TCCTCCGCTTATTGATATGC"

r1_its1 = []
r1_its2 = []
r2_its1 = []
r2_its2 = []
undetermined = []

readiter = SeqIO.parse(open(inpath), "fastq")

for rec1, rec2 in izip_longest(readiter, readiter):
	if primer_its1f in rec1.seq:
		r1_its1.append(rec1)
	if primer_its2r in rec1.seq:
		r2_its1.append(rec1)
	if primer_its3f in rec1.seq:
		r1_its1.append(rec1)
	if primer_its4r in rec1.seq:
		r2_its2.append(rec1)
	else:
		undetermined.append(rec1, rec2)
		continue

        if primer_its1f in rec2.seq:
                r1_its1.append(rec1)
        if primer_its2r in rec2.seq:
                r2_its1.append(rec1)
        if primer_its3f in rec2.seq:
                r1_its1.append(rec1)
        if primer_its4r in rec2.seq:
                r2_its2.append(rec1)

outITS1_R1 = open("ITS1_R1.fastq")
outITS1_R2 = open("ITS1_R2.fastq")
outITS2_R1 = open("ITS2_R1.fastq")
outITS2_R2 = open("ITS2_R2.fastq")

SeqIO.write(r1_its1, outITS1_R1, "fastq")
outITS1_R1.close()

SeqIO.write(r2_its1, outITS1_R2, "fastq")
outITS1_R2.close()

SeqIO.write(r1_its2, outITS2_R1, "fastq")
outITS2_R1.close()

SeqIO.write(r2_its2, outITS2_R2, "fastq")
outITS2_R2.close()

