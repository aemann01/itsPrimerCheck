#!/usr/bin/env python
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys

'''Pulls ITS regions from genbank files. Useage: python genbankSeqExtract.py in.gb > out.fa'''

gb_file = sys.argv[1]

for gb_record in SeqIO.parse(open(gb_file, "r"), "genbank"):
	for feature in gb_record.features:
		if feature.type == "misc_RNA" or feature.type == "RNA":
			if 'product' in feature.qualifiers: 
				if 'internal transcribed' in feature.qualifiers['product'][0]:
					start = feature.location.start
					end = feature.location.end
					print(">%s|%s|%s|%i\n%s" % (gb_record.description, gb_record.name, feature.qualifiers['product'][0], len(gb_record.seq[start:end]), gb_record.seq[start:end])) 
			elif 'locus_tag' in feature.qualifiers:
				if 'ITS' in feature.qualifiers['locus_tag'][0]:
					start = feature.location.start
					end = feature.location.end
					print(">%s|%s|%s|%i\n%s" % (gb_record.description, gb_record.name, feature.qualifiers['locus_tag'][0], len(gb_record.seq[start:end]), gb_record.seq[start:end]))
			elif 'gene' in feature.qualifiers:
                                if 'ITS' in feature.qualifiers['gene'][0]:
                                        start = feature.location.start
                                        end = feature.location.end
                                        print(">%s|%s|%s|%i\n%s" % (gb_record.description, gb_record.name, feature.qualifiers['gene'][0], len(gb_record.seq[start:end]), gb_record.seq[start:end]))
			
			#this doesn't work at the moment, will only print if the called feature qualifiers are not there
			#else:
			#	print("Unable to file ITS entry for %s, filename %s" % (feature.organism, gb_file))
