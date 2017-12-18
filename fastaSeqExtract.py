#!/usr/bin/env python

'''For filtering out RNA seqs out of fasta files. Useage: python fastaSeqExtract.py in.fasta out.fasta'''

import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

rnafa = []

for record in SeqIO.parse(sys.argv[1], "fasta"):
	if "S ribosomal RNA " and not "mRNA" in record.description:
		rnafa.append(record)
	else:
		continue

with open(sys.argv[2], "w") as handle:
	SeqIO.write(rnafa, handle, "fasta")
