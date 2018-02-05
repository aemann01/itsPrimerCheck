#Not working version -- choppy from ipython

import sys
import pandas as pd
from Bio.Blast import NCBIXML

for rec in NCBIXML.parse(open("its1.blast.out.xml")):
	for alignment in rec.alignments:
		for hsp in alignment.hsps:
			print("%s\t%i\t%i\t%f\t%s" % (rec.query, hsp.identities, hsp.align_length, (hsp.identities/hsp.align_length), alignment.title))

#merge length and blastout
#tax = pd.read_csv("genera.its2.txt", sep="\t")
#len = pd.read_csv("../rep_set/its2rep97.eukOnly.lens", sep="\t")
#merge = pd.merge(tax, len)
#merge['PID'] = merge

tax1 = pd.read_csv("genera.its1.txt", sep="\t")
tax2 = pd.read_csv("genera.its2.txt", sep="\t")

len1 = pd.read_csv("../rep_set/its1rep97.eukOnly.lens", sep="\t")
len2 = pd.read_csv("../rep_set/its2rep97.eukOnly.lens", sep="\t")

merge1 = pd.merge(tax1, len1)
merge2 = pd.merge(tax2, len2)

merge1['PID'] = merge1['IDENT']/merge1['LENGTH']
merge2['PID'] = merge2['IDENT']/merge2['LENGTH']

#now get top per genera and order the output
merge1top = merge1.groupby(['OTU', 'TAX'])['PID'].max().reset_index().sort_values(['OTU', 'PID'], ascending=False).set_index(['OTU', 'TAX']).reset_index()
merge2top = merge2.groupby(['OTU', 'TAX'])['PID'].max().reset_index().sort_values(['OTU', 'PID'], ascending=False).set_index(['OTU', 'TAX']).reset_index()

shortlist1 = merge1top.set_index(['TAX']).groupby('OTU')['PID'].nlargest(5).reset_index()
shortlist2 = merge2top.set_index(['TAX']).groupby('OTU')['PID'].nlargest(5).reset_index()

#now get the rest of the metadata back merged
pd.merge(shortlist1, merge1, how='left').drop_duplicates(shortlist1.columns.difference(['ACCESS'])).to_csv("generaShortList.its1.txt", sep="\t", index=False)
pd.merge(shortlist2, merge2, how='left').drop_duplicates(shortlist2.columns.difference(['ACCESS'])).to_csv("generaShortList.its2.txt", sep="\t", index=False)






