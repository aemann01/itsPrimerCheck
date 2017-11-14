#!/usr/bin/env python
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys

'''Pulls ITS regions from genbank files. Useage: python genbankSeqExtract.py in.gb > out.fa'''

gb_file = sys.argv[1]

for gb_record in SeqIO.parse(open(gb_file, "r"), "genbank"):
	for feature in gb_record.features:
		if(feature.type == "misc_RNA"):
			if 'ITS' in feature.qualifiers['locus_tag'][0]:
				start = feature.location.start
				end = feature.location.end
				print(">%s|%s|%s|%i\n%s") % (gb_record.description, gb_record.name, feature.qualifiers['locus_tag'][0], len(gb_record.seq[start:end]), gb_record.seq[start:end])			
