#!/usr/bin/env python

#by Bianca De Sanctis <bdd28@cam.ac.uk>

import sys
from Bio import AlignIO
from Bio.Align import AlignInfo

alignment = AlignIO.read("extended_mammoth_published_dna.afa", "fasta")
summary_align = AlignInfo.SummaryInfo(alignment)
consensus = summary_align.dumb_consensus()
print(len(consensus))

with open('extended_mammoth_consensus.fa', 'a') as the_file:
	the_file.write(str(consensus))

