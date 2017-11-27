#!/usr/bin/env/python3

import sys
from pandas import DataFrame
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

gb_file = sys.argv[1]

ITS1_fa = []
ITS2_fa = []
5S_fa = []
stat_tb = []

for record in SeqIO.parse(open(gb_file, "r"), "genbank"):
	taxString = record.//////
	speciesID = record.//////

	#gather coordinate information
	if feature.type == "RNA" and not "mRNA":
		if "product" in feature.qualifiers:
			if "18S ribosomal" in feature.qualifiers['product'][0]:
				if feature.location.end > feature.location.start:
					18S_end = feature.location.end
				else:
					18S_end = feature.location.start

			if "5S ribosomal" or "5.8S ribosomal" in feature.qualifiers['product'][0]:
				if feature.location.end > feature.location.start:
					5S_end = feature.location.end
					5S_start = feature.location.start
				else:
					5S_end = feature.location.start
					5S_start = feature.location.end

			if "28S ribosomal" or "26S ribosomal" in feature.qualifiers['product'][0]:
				if feature.location.start < feature.location.end:
					26S_start = feature.location.start
				else:
					26S_start = feature.location.end

			ITS1_fa.append

		ITS1_len = 5S_start-18S_end
		ITS2_len = 26S_start-5S_end
		5S_len = 5S_start-5S_end

#population statistics table
stat_tb.append([speciesID, taxString, 18S_end, 5S_start, 5S_end, 26S_start, ITS1_len, ITS2_len, 5S_len])
df = DataFrame(stat_tb)
df.columns = ["speciesID", "taxString", "18S_end", "5S_start", "5S_end", "26S_start", "ITS1_len", "ITS2_len", "5S_len"]


