import sys
import gzip
import pandas as pd
import glob
import uuid
from pandas import DataFrame
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

files = glob.glob(".gz")

es_query = ["18S rRNA", "18S ribosomal"]
fs_query = ["5S rRNA", "5S ribosomal", "5.8S rRNA", "5.8S ribosomal"]
ts_query = ["26S rRNA", "26S ribosomal", "28S rRNA", "28S ribosomal"]
stat_tb = []

for f in files:
	for record in SeqIO.parse(gzip.open(f, "r"), "genbank"):
		for feature in record.features:
			if feature.type == "RNA" or "CDS" and "product" in feature.qualifiers:
				if any(x in feature.qualifiers["product"][0] for x in es_query):
					if feature.location.end > feature.location.start:
						Es_end = feature.location.end
					elif feature.location.end < feature.location.start:
						Es_end = feature.location.start
				if any(x in feature.qualifiers["product"][0] for x in fs_query):
					if feature.location.end > feature.location.start:
						Fs_end = feature.location.end
						Fs_start = feature.location.start
					elif feature.location.end < feature.location.start:
						Fs_end = feature.location.start
						Fs_start = feature.location.end
				if any(x in feature.qualifiers["product"][0] for x in ts_query):
					if feature.location.end > feature.location.start:
						Ts_start = feature.location.start
					elif feature.location.end < feature.location.start:
						Ts_start = feature.location.end
				try:
					Es_end
				except NameError:
					Es_end = "NA"
				try:
					Fs_end
				except NameError:
					Fs_end = "NA"
				try:
					Fs_start
				except NameError:
					Fs_start = "NA"
				try:
					Ts_end
				except NameError:
					Ts_end = "NA"


	stat_tb.append([f, Es_end, Fs_start, Fs_end, Ts_start])
	print(stat_tb)
df = DataFrame(stat_tb)

outfile = str(uuid.uuid4())

with open(outfile, "w") as f:
	for line in df:
		f.write(line)

