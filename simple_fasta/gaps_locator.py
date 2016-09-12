#!/usr/bin/env python
### Counts N (actual N, not A,T,G,C)
### prints out placement of N strings
### for each seq in fasta file
import string
import sys
import Bio
import re
from Bio import SeqIO
from StringIO import StringIO

fastafile = sys.argv[1]
pattern = re.compile(r'([N]+)') 
handle = open(fastafile)
for seq_record in SeqIO.parse(handle, "fasta"):
	matches = re.finditer(pattern, str(seq_record.seq))
	Ncnt = seq_record.seq.count("N")
	print str(seq_record.id) +' total N count= '+ str(Ncnt)
	if Ncnt > 1:
		print 'gaps at:'
		print [str(m.start(1))+'-'+str(m.end(1)) for m in matches]	
handle.close()


