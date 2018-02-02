import sys
from Bio.Blast import NCBIXML

for rec in NCBIXML.parse(open("its2.blast.out.xml")):
	for alignment in rec.alignments:
		for hsp in alignment.hsps:
			print("%s\t%s\t%s\t%i\f%i\t%f" % (rec.query, alignment.title, alignment.length, hsp.identities, alignment.length, (hsp.identities/alignment.length)))
