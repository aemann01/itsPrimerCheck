import sys
import gzip
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

data = []

for record in SeqIO.parse(open("repSet_blastResults.gb", "r"), "genbank"):
     for feature in record.features:
         accession_num = record.name
         taxa_string = ";".join(record.annotations["taxonomy"])
         source_org = record.annotations["source"].replace(" ", "_")
         data.append([accession_num, source_org, taxa_string])
     df = pd.DataFrame(data)
     df.columns = ['accession', 'source', 'taxonomy']
df.to_csv("taxonomy.txt", sep="\t", index=False)

