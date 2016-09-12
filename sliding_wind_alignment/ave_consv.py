#!/usr/bin/env python

'''
DESCRIPTION: 
retrieves nt conservation from 0 - 1 (1 being perfectly conserved)
at each site from a pair of aligned sequences.

USAGE: 
python scriptname.py infile.fasta outfile.tab windowsize
where infile == pairwise alignment, FASTA format, only two sequences
and windowsize is number only (eg 100, 200, 500, 1000)

NOTES: 
Treats gaps (- and N) as 0, thus not conserved 
'''

import sys
import os
from Bio import AlignIO
from collections import Counter

def get_per_cons(AlignFile):
    perIDlist = []
    Align = AlignIO.read(AlignFile, "fasta")
    Sites = len(Align[0])
    Denom = float(len(Align))
    entriesToRemove = ('n', 'N', '-')

    for x in range(Sites):
        SiteClctn = Counter(Align[:,x])

        for k in entriesToRemove:
            del SiteClctn[k]

        ComChar = SiteClctn.most_common(1)
        if len(ComChar) == 0:
            perIDlist.append(0)
        if len(ComChar) >= 1:
            if int(ComChar[0][1]) == 2:
                perIDlist.append(1)
            if int(ComChar[0][1]) == 1:
                perIDlist.append(0)
            if len(ComChar) > 2 or int(ComChar[0][1]) < 0:
                print 'ERROR, check Aligned file and output'
                sys.exit('Exiting script...')
    return (perIDlist)

def slidingWindow(perIDlist, WinSize):
    #this function written by user xguse, see http://en.gravatar.com/xguse
    Step=1
    WinVals = []

    try: 
        it = iter(perIDlist)
    except TypeError:
        raise Exception("**ERROR** sequence must be iterable.")
    if not ((type(WinSize) == type(0)) and (type(Step) == type(0))):
        raise Exception("**ERROR** type(winSize) and type(step) must be int.")
    if Step > WinSize:
        raise Exception("**ERROR** step must not be larger than winSize.")
    if WinSize > len(perIDlist):
        raise Exception('''**ERROR** winSize must not 
                        be larger than sequence length.''')
 
    #Pre-compute number of chunks to emit
    numOfChunks = ((len(perIDlist)-WinSize)/Step)+1
 
    #Chunking
    for i in range(0,numOfChunks*Step,Step):
        WinVals.append(reduce(lambda x, y: x + y, perIDlist[i:i+WinSize]) 
                        / float(WinSize))
    return (WinVals)

def print_out(WinVals, OutFile):
    OutHandle = open(str(OutFile), 'w')
    for RowNum, WinVal in enumerate(WinVals, start=1):
        OutHandle.write('%s\t%s\n' % (RowNum, WinVal))
    OutHandle.close()

InFile = sys.argv[1]
OutFile = sys.argv[2]
WinSize = sys.argv[3]

perIDlist = get_per_cons(InFile)
WinVals = slidingWindow(perIDlist, int(WinSize))
print_out(WinVals, OutFile)
